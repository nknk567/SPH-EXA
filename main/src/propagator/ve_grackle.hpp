

#pragma once

#include <variant>

#include "cstone/fields/field_get.hpp"
#include "sph/particles_data.hpp"
#include "sph/sph.hpp"

#include "ipropagator.hpp"
#include "gravity_wrapper.hpp"

#include "cooling/cooler.hpp"
#include "cooling/eos_cooling.hpp"

namespace sphexa
{

using namespace sph;
using util::FieldList;

template<bool avClean, class DomainType, class DataType>
class VeCooling : public Propagator<DomainType, DataType>
{
protected:
    using Base = Propagator<DomainType, DataType>;
    using Base::timer;

    using T             = typename DataType::RealType;
    using KeyType       = typename DataType::KeyType;
    using Tmass         = typename DataType::HydroData::Tmass;
    using ChemData      = typename DataType::ChemData;
    using MultipoleType = ryoanji::CartesianQuadrupole<Tmass>;

    using Acc       = typename DataType::AcceleratorType;
    using MHolder_t = typename cstone::AccelSwitchType<Acc, MultipoleHolderCpu, MultipoleHolderGpu>::template type<
        MultipoleType, DomainType, typename DataType::HydroData>;

    MHolder_t          mHolder_;
    cooling::Cooler<T> cooling_data;

    /*! @brief the list of conserved particles fields with values preserved between iterations
     *
     * x, y, z, h and m are automatically considered conserved and must not be specified in this list
     */
    using ConservedFields = FieldList<"u", "vx", "vy", "vz", "x_m1", "y_m1", "z_m1", "du_m1", "alpha">;

    //! @brief list of dependent fields, these may be used as scratch space during domain sync
    using DependentFields_ = FieldList<"prho", "c", "ax", "ay", "az", "du", "c11", "c12", "c13", "c22", "c23", "c33",
                                       "xm", "kx", "nc", "rho", "p">;

    //! @brief velocity gradient fields will only be allocated when avClean is true
    using GradVFields = FieldList<"dV11", "dV12", "dV13", "dV22", "dV23", "dV33">;

    //! @brief what will be allocated based AV cleaning choice
    using DependentFields =
        std::conditional_t<avClean, decltype(DependentFields_{} + GradVFields{}), decltype(DependentFields_{})>;
    using CoolingFields = typename util::MakeFieldList<ChemData>::Fields;

public:
    VeCooling(std::ostream& output, size_t rank)
        : Base(output, rank)
    {
        if (avClean && rank == 0) { std::cout << "AV cleaning is activated" << std::endl; }
    }

    std::vector<std::string> conservedFields() const override
    {
        std::vector<std::string> ret{"x", "y", "z", "h", "m"};
        for_each_tuple([&ret](auto f) { ret.push_back(f.value); }, make_tuple(ConservedFields{}));
        return ret;
    }

    void load(const std::string& initCond, MPI_Comm comm) override
    {
        if (initCond == "evrard-cooling")
        {
            BuiltinWriter attributeSetter(evrardCoolingConstants());
            cooling_data.loadOrStoreAttributes(&attributeSetter);
        }
        else
        {
            std::string path = strBeforeSign(initCond, ":");
            if (std::filesystem::exists(path))
            {
                std::unique_ptr<IFileReader> reader;
                reader = std::make_unique<H5PartReader>(comm);
                reader->setStep(path, -1);
                cooling_data.loadOrStoreAttributes(reader.get());
                reader->closeStep();
            }
            else
                throw std::runtime_error("");
        }
        cooling_data.init(0);
    }

    void save(IFileWriter* writer) override { cooling_data.loadOrStoreAttributes(writer); }

    void activateFields(DataType& simData) override
    {
        auto& d = simData.hydro;
        //! @brief Fields accessed in domain sync (x,y,z,h,m,keys) are not part of extensible lists.
        d.setConserved("x", "y", "z", "h", "m");
        d.setDependent("keys");
        std::apply([&d](auto... f) { d.setConserved(f.value...); }, make_tuple(ConservedFields{}));
        std::apply([&d](auto... f) { d.setDependent(f.value...); }, make_tuple(DependentFields{}));

        d.devData.setConserved("x", "y", "z", "h", "m");
        d.devData.setDependent("keys");
        std::apply([&d](auto... f) { d.devData.setConserved(f.value...); }, make_tuple(ConservedFields{}));
        std::apply([&d](auto... f) { d.devData.setDependent(f.value...); }, make_tuple(DependentFields{}));
    }

    void sync(DomainType& domain, DataType& simData) override
    {
        auto& d = simData.hydro;
        if (d.g != 0.0)
        {
            domain.syncGrav(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                            std::tuple_cat(get<ConservedFields>(d), get<CoolingFields>(simData.chem)),
                            get<DependentFields>(d));
        }
        else
        {
            domain.sync(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                            std::tuple_cat(get<ConservedFields>(d), get<CoolingFields>(simData.chem)),
                            get<DependentFields>(d));
        }
        d.treeView = domain.octreeProperties();
    }

    void computeForces(DomainType& domain, DataType& simData)
    {
        timer.start();
        sync(domain, simData);
        timer.step("domain::sync");

        auto& d = simData.hydro;
        d.resize(domain.nParticlesWithHalos());
        resizeNeighbors(d, domain.nParticles() * d.ngmax);
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        transferToHost(d, first, first + 1, {"m"});
        fill(get<"m">(d), 0, first, d.m[first]);
        fill(get<"m">(d), last, domain.nParticlesWithHalos(), d.m[first]);

        findNeighborsSfc(first, last, d, domain.box());
        timer.step("FindNeighbors");

        computeXMass(first, last, d, domain.box());
        timer.step("XMass");
        domain.exchangeHalos(std::tie(get<"xm">(d)), get<"ax">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        d.release("ax");
        d.devData.release("ax");
        d.acquire("gradh");
        d.devData.acquire("gradh");
        computeVeDefGradh(first, last, d, domain.box());
        timer.step("Normalization & Gradh");

        // computeEOS(first, last, d);
        eos_cooling_ve(first, last, d, simData.chem, cooling_data);
        timer.step("EquationOfState");

        domain.exchangeHalos(get<"vx", "vy", "vz", "prho", "c", "kx">(d), get<"gradh">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        d.release("gradh", "ay");
        d.devData.release("gradh", "ay");
        d.acquire("divv", "curlv");
        d.devData.acquire("divv", "curlv");
        computeIadDivvCurlv(first, last, d, domain.box());
        d.minDtRho = rhoTimestep(first, last, d);
        timer.step("IadVelocityDivCurl");

        domain.exchangeHalos(get<"c11", "c12", "c13", "c22", "c23", "c33", "divv">(d), get<"az">(d), get<"du">(d));
        timer.step("mpi::synchronizeHalos");

        computeAVswitches(first, last, d, domain.box());
        timer.step("AVswitches");

        if (avClean)
        {
            domain.exchangeHalos(get<"dV11", "dV12", "dV22", "dV23", "dV33", "alpha">(d), get<"az">(d), get<"du">(d));
        }
        else { domain.exchangeHalos(std::tie(get<"alpha">(d)), get<"az">(d), get<"du">(d)); }
        timer.step("mpi::synchronizeHalos");

        d.release("divv", "curlv");
        d.devData.release("divv", "curlv");
        d.acquire("ax", "ay");
        d.devData.acquire("ax", "ay");
        computeMomentumEnergy<avClean>(first, last, d, domain.box());
        timer.step("MomentumAndEnergy");

        if (d.g != 0.0)
        {
            mHolder_.upsweep(d, domain);
            timer.step("Upsweep");
            mHolder_.traverse(d, domain);
            timer.step("Gravity");
        }
    }

    void step(DomainType& domain, DataType& simData) override
    {
        computeForces(domain, simData);

        auto&  d     = simData.hydro;
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        auto minDtCooling = cooling::coolingTimestep(first, last, d, cooling_data, simData.chem);
        computeTimestep(first, last, d, minDtCooling);
        timer.step("Timestep");

#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; i++)
        {
            T u_old  = d.u[i];
            T u_cool = d.u[i];
            cooling_data.cool_particle(d.minDt, d.rho[i], u_cool,
                                       cstone::getPointers(get<CoolingFields>(simData.chem), i));
            const T du = (u_cool - u_old) / d.minDt;
            d.du[i] += du;
        }
        timer.step("GRACKLE chemistry and cooling");

        computePositions(first, last, d, domain.box());
        timer.step("UpdateQuantities");
        updateSmoothingLength(first, last, d);
        timer.step("UpdateSmoothingLength");

        timer.stop();
    }

    void saveFields(IFileWriter* writer, size_t first, size_t last, DataType& simData,
                    const cstone::Box<T>& box) override
    {
        auto&            d             = simData.hydro;
        auto             fieldPointers = d.data();
        std::vector<int> outputFields  = d.outputFieldIndices;

        auto output = [&]()
        {
            for (int i = int(outputFields.size()) - 1; i >= 0; --i)
            {
                int fidx = outputFields[i];
                if (d.isAllocated(fidx))
                {
                    int column = std::find(d.outputFieldIndices.begin(), d.outputFieldIndices.end(), fidx) -
                                 d.outputFieldIndices.begin();
                    transferToHost(d, first, last, {d.fieldNames[fidx]});
                    std::visit([writer, c = column, key = d.fieldNames[fidx]](auto field)
                               { writer->writeField(key, field->data(), c); },
                               fieldPointers[fidx]);
                    outputFields.erase(outputFields.begin() + i);
                }
            }
        };

        // first output pass: write everything allocated at the end of the step
        output();

        d.release("ax", "ay", "az");
        d.devData.release("ax", "ay", "az");

        // second output pass: write temporary quantities produced by the EOS
        d.acquire("rho", "p", "gradh");
        d.devData.acquire("rho", "p", "gradh");
        computeEOS(first, last, d);
        output();
        d.devData.release("rho", "p", "gradh");
        d.release("rho", "p", "gradh");

        // third output pass: curlv and divv
        d.acquire("divv", "curlv");
        d.devData.acquire("divv", "curlv");
        if (!outputFields.empty()) { computeIadDivvCurlv(first, last, d, box); }
        output();
        d.release("divv", "curlv");
        d.devData.release("divv", "curlv");

        d.acquire("ax", "ay", "az");
        d.devData.acquire("ax", "ay", "az");

        if (!outputFields.empty() && Base::rank_ == 0)
        {
            std::cout << "WARNING: the following fields are not in use and therefore not output: ";
            for (int fidx = 0; fidx < outputFields.size() - 1; ++fidx)
            {
                std::cout << d.fieldNames[fidx] << ",";
            }
            std::cout << d.fieldNames[outputFields.back()] << std::endl;
        }
    }
};

} // namespace sphexa
