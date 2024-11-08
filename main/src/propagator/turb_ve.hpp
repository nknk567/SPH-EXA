/*
 * MIT License
 *
 * Copyright (c) 2021 CSCS, ETH Zurich
 *               2021 University of Basel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*! @file
 * @brief A Propagator class for modern SPH with generalized volume elements
 *
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

#pragma once

#include <filesystem>
#include <sstream>
#include <variant>

#include "cstone/util/constexpr_string.hpp"
#include "cstone/fields/field_get.hpp"
#include "sph/sph.hpp"
#include "sph/hydro_turb/turbulence_data.hpp"

#include "ve_hydro.hpp"
#include "gravity_wrapper.hpp"

#include "adaptive_timestep_helper.hpp"

namespace sphexa
{

using namespace sph;
//! @brief VE hydro propagator that adds turbulence stirring to the acceleration prior to position update
template<bool avClean, class DomainType, class DataType>
class TurbVeBdtProp final : public Propagator<DomainType, DataType>
{
protected:
    using Base = Propagator<DomainType, DataType>;
    using Base::pmReader;
    using Base::timer;

    using T             = typename DataType::RealType;
    using KeyType       = typename DataType::KeyType;
    using Tmass         = typename DataType::HydroData::Tmass;
    using MultipoleType = ryoanji::CartesianQuadrupole<Tmass>;

    using Acc       = typename DataType::AcceleratorType;
    using MHolder_t = typename cstone::AccelSwitchType<Acc, MultipoleHolderCpu, MultipoleHolderGpu>::template type<
        MultipoleType, DomainType, typename DataType::HydroData>;

    AdaptiveTimestepHelper<DomainType, DataType> timestepHelper;
    MHolder_t                                    mHolder_;

    sph::TurbulenceData<typename DataType::RealType, typename DataType::AcceleratorType> turbulenceData;
    using ConservedFields = FieldList<"c", "vx", "vy", "vz", "x_m1", "y_m1", "z_m1", "du_m1", "alpha", "rung">;

    //! @brief list of dependent fields, these may be used as scratch space during domain sync
    using DependentFields_ = FieldList<"ax", "ay", "az", "prho", "du", "c11", "c12", "c13", "c22", "c23", "c33", "xm",
                                       "kx", "nc", "divv", "gradh", "u">;

    //! @brief velocity gradient fields will only be allocated when avClean is true
    using GradVFields = FieldList<"dV11", "dV12", "dV13", "dV22", "dV23", "dV33">;

    //! @brief what will be allocated based AV cleaning choice
    using DependentFields =
        std::conditional_t<avClean, decltype(DependentFields_{} + GradVFields{}), decltype(DependentFields_{})>;

public:
    TurbVeBdtProp(std::ostream& output, size_t rank, const InitSettings& settings)
        : Base(output, rank)
        , timestepHelper(settings)
        , turbulenceData(settings, rank == 0)
    {
    }

    std::vector<std::string> conservedFields() const override
    {
        std::vector<std::string> ret{"x", "y", "z", "h", "m"};
        for_each_tuple([&ret](auto f) { ret.push_back(f.value); }, make_tuple(ConservedFields{}));
        return ret;
    }

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

    void save(IFileWriter* writer) override
    {
        timestepHelper.save(writer);
        Base::save(writer);
        turbulenceData.loadOrStore(writer);
    }

    void load(const std::string& initCond, IFileReader* reader) override
    {
        timestepHelper.load(initCond, reader);
        Base::load(initCond, reader);

        int         step = numberAfterSign(initCond, ":");
        std::string path = removeModifiers(initCond);
        // The file does not exist, we're starting from scratch. Nothing to do.
        if (!std::filesystem::exists(path)) { return; }

        reader->setStep(path, step, FileMode::independent);
        turbulenceData.loadOrStore(reader);

        if (Base::rank_ == 0) { std::cout << "Restored turbulence state from " << path << ":" << step << std::endl; }
        reader->closeStep();
    }

    void sync(DomainType& domain, DataType& simData) override
    {
        timestepHelper.template syncAdaptiveDt<ConservedFields, DependentFields>(domain, simData);
    }

    void computeForces(DomainType& domain, DataType& simData) override
    {
        timer.start();
        pmReader.start();
        sync(domain, simData);
        timer.step("domain::sync");

        auto&  d     = simData.hydro;
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        transferToHost(d, first, first + 1, {"m"});
        fill(get<"m">(d), 0, first, d.m[first]);
        fill(get<"m">(d), last, domain.nParticlesWithHalos(), d.m[first]);

        findNeighborsSfc(first, last, d, domain.box());
        timer.step("FindNeighbors");
        pmReader.step();

        computeXMass(timestepHelper.activeRungs_, d, domain.box());
        timer.step("XMass");
        domain.exchangeHalos(std::tie(get<"xm">(d)), get<"keys">(d), timestepHelper.haloRecvScratch);
        timer.step("mpi::synchronizeHalos");

        computeVeDefGradh(timestepHelper.activeRungs_, d, domain.box());
        timer.step("Normalization & Gradh");

        computeIsothermalEOS(first, last, d);
        timer.step("EquationOfState");

        domain.exchangeHalos(get<"vx", "vy", "vz", "prho", "c", "kx">(d), get<"keys">(d),
                             timestepHelper.haloRecvScratch);
        timer.step("mpi::synchronizeHalos");

        computeIadDivvCurlv(timestepHelper.activeRungs_, d, domain.box());
        groupDivvTimestep(timestepHelper.activeRungs_, rawPtr(timestepHelper.groupDt_), d);
        timer.step("IadVelocityDivCurl");

        domain.exchangeHalos(get<"c11", "c12", "c13", "c22", "c23", "c33", "divv">(d), get<"keys">(d),
                             timestepHelper.haloRecvScratch);
        timer.step("mpi::synchronizeHalos");

        computeAVswitches(timestepHelper.activeRungs_, d, domain.box());
        timer.step("AVswitches");

        if (avClean)
        {
            domain.exchangeHalos(get<"dV11", "dV12", "dV22", "dV23", "dV33", "alpha">(d), get<"keys">(d),
                                 timestepHelper.haloRecvScratch);
        }
        else { domain.exchangeHalos(std::tie(get<"alpha">(d)), get<"keys">(d), timestepHelper.haloRecvScratch); }
        timer.step("mpi::synchronizeHalos");

        computeMomentumEnergy<avClean>(timestepHelper.activeRungs_, rawPtr(timestepHelper.groupDt_), d, domain.box());
        timer.step("MomentumAndEnergy");
        pmReader.step();

        if (d.g != 0.0)
        {
            //            bool      isNewHierarchy = timestepHelper.activeRung(timestepHelper.timestep_.substep,
            //            timestepHelper.timestep_.numRungs) == 0; GroupView gravGroup      = isNewHierarchy ?
            //            mHolder_.computeSpatialGroups(d, domain) : timestepHelper.activeRungs_;
            GroupView gravGroup = timestepHelper.getGravGroupView(mHolder_, domain, simData);

            mHolder_.upsweep(d, domain);
            timer.step("Upsweep");
            pmReader.step();
            mHolder_.traverse(gravGroup, d, domain);
            timer.step("Gravity");
            pmReader.step();
        }
        groupAccTimestep(timestepHelper.activeRungs_, rawPtr(timestepHelper.groupDt_), d);
        driveTurbulence(timestepHelper.activeRungs_, simData.hydro, turbulenceData);
        timer.step("Turbulence Stirring");
    }

    void integrate(DomainType& domain, DataType& simData) override
    {
        timestepHelper.integrate(domain, simData, timer, Base::rank_);
    }
};

template<bool avClean, class DomainType, class DataType>
class TurbVeProp final : public HydroVeProp<avClean, DomainType, DataType>
{
    using Base = HydroVeProp<avClean, DomainType, DataType>;
    using Base::rank_;
    using Base::timer;

    sph::TurbulenceData<typename DataType::RealType, typename DataType::AcceleratorType> turbulenceData;

public:
    TurbVeProp(std::ostream& output, size_t rank, const InitSettings& settings)
        : Base(output, rank)
        , turbulenceData(settings, rank == 0)
    {
    }

    void computeForces(DomainType& domain, DataType& simData) override
    {
        Base::computeForces(domain, simData);
        driveTurbulence(Base::groups_.view(), simData.hydro, turbulenceData);
        timer.step("Turbulence Stirring");
    }

    void save(IFileWriter* writer) override { turbulenceData.loadOrStore(writer); }

    void load(const std::string& initCond, IFileReader* reader) override
    {
        int         step = numberAfterSign(initCond, ":");
        std::string path = removeModifiers(initCond);
        // The file does not exist, we're starting from scratch. Nothing to do.
        if (!std::filesystem::exists(path)) { return; }

        reader->setStep(path, step, FileMode::independent);
        turbulenceData.loadOrStore(reader);

        if (rank_ == 0) { std::cout << "Restored turbulence state from " << path << ":" << step << std::endl; }
        reader->closeStep();
    }
};

} // namespace sphexa
