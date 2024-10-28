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

namespace sphexa
{

using namespace sph;

//! @brief VE hydro propagator that adds turbulence stirring to the acceleration prior to position update
template<bool avClean, class DomainType, class DataType>
class TurbVeBdtProp final : public HydroVeBdtProp<avClean, DomainType, DataType>
{
    using Base = HydroVeBdtProp<avClean, DomainType, DataType>;
    using Base::rank_;
    using Base::timer;

    sph::TurbulenceData<typename DataType::RealType, typename DataType::AcceleratorType> turbulenceData;

public:
    TurbVeBdtProp(std::ostream& output, size_t rank, const InitSettings& settings)
        : Base(output, rank, settings)
        , turbulenceData(settings, rank == 0)
    {
    }

    void computeForces(DomainType& domain, DataType& simData) override
    {
        Base::computeForces(domain, simData);
        driveTurbulence(Base::activeRungs_, simData.hydro, turbulenceData);
        timer.step("Turbulence Stirring");
    }

    void save(IFileWriter* writer) override
    {
        Base::save(writer);
        turbulenceData.loadOrStore(writer);
    }

    void load(const std::string& initCond, IFileReader* reader) override
    {
        Base::load(initCond, reader);

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

template<bool avClean, class DomainType, class DataType>
class TurbVeProp final : public HydroVeProp<avClean, DomainType, DataType>
{
    using Base = HydroVeProp<avClean, DomainType, DataType>;
    using Base::rank_;
    using Base::timer;
    using Base::pmReader;
    using Base::mHolder_;

    sph::TurbulenceData<typename DataType::RealType, typename DataType::AcceleratorType> turbulenceData;
    using ConservedFields = FieldList<"temp", "vx", "vy", "vz", "x_m1", "y_m1", "z_m1", "du_m1", "alpha", "c">;

    //! @brief list of dependent fields, these may be used as scratch space during domain sync
    using DependentFields_ =
        FieldList<"ax", "ay", "az", "prho", "du", "c11", "c12", "c13", "c22", "c23", "c33", "xm", "kx", "nc">;

    //! @brief velocity gradient fields will only be allocated when avClean is true
    using GradVFields = FieldList<"dV11", "dV12", "dV13", "dV22", "dV23", "dV33">;

    //! @brief what will be allocated based AV cleaning choice
    using DependentFields =
        std::conditional_t<avClean, decltype(DependentFields_{} + GradVFields{}), decltype(DependentFields_{})>;

public:
    TurbVeProp(std::ostream& output, size_t rank, const InitSettings& settings)
        : Base(output, rank)
        , turbulenceData(settings, rank == 0)
    {
    }
    void sync(DomainType& domain, DataType& simData) override
    {
        auto& d = simData.hydro;
        if (d.g != 0.0)
        {
            domain.syncGrav(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                            get<ConservedFields>(d), get<DependentFields>(d));
        }
        else
        {
            domain.sync(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d),
                        std::tuple_cat(std::tie(get<"m">(d)), get<ConservedFields>(d)), get<DependentFields>(d));
        }
        d.treeView = domain.octreeProperties();
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


    void computeForces(DomainType& domain, DataType& simData) override
    {
        timer.start();
        pmReader.start();
        sync(domain, simData);
        timer.step("domain::sync");

        auto& d = simData.hydro;
        d.resizeAcc(domain.nParticlesWithHalos());
        resizeNeighbors(d, domain.nParticles() * d.ngmax);
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        transferToHost(d, first, first + 1, {"m"});
        fill(get<"m">(d), 0, first, d.m[first]);
        fill(get<"m">(d), last, domain.nParticlesWithHalos(), d.m[first]);

        findNeighborsSfc(first, last, d, domain.box());
        computeGroups(first, last, d, domain.box(), Base::groups_);
        timer.step("FindNeighbors");
        pmReader.step();

        computeXMass(Base::groups_.view(), d, domain.box());
        timer.step("XMass");
        domain.exchangeHalos(std::tie(get<"xm">(d)), get<"ax">(d), get<"keys">(d));
        timer.step("mpi::synchronizeHalos");

        release(d, "ay");
        acquire(d, "gradh");
        computeVeDefGradh(Base::groups_.view(), d, domain.box());
        timer.step("Normalization & Gradh");

        computeIsothermalEOS(first, last, d);
        timer.step("EquationOfState");

        domain.exchangeHalos(get<"vx", "vy", "vz", "prho", "c", "kx">(d), get<"ax">(d), get<"keys">(d));
        timer.step("mpi::synchronizeHalos");

        release(d, "gradh", "az");
        acquire(d, "divv", "curlv");
        computeIadDivvCurlv(Base::groups_.view(), d, domain.box());
        d.minDtRho = rhoTimestep(first, last, d);
        timer.step("IadVelocityDivCurl");

        domain.exchangeHalos(get<"c11", "c12", "c13", "c22", "c23", "c33", "divv">(d), get<"ax">(d), get<"keys">(d));
        timer.step("mpi::synchronizeHalos");

        computeAVswitches(Base::groups_.view(), d, domain.box());
        timer.step("AVswitches");

        if (avClean)
        {
            domain.exchangeHalos(get<"dV11", "dV12", "dV22", "dV23", "dV33", "alpha">(d), get<"ax">(d), get<"keys">(d));
        }
        else { domain.exchangeHalos(std::tie(get<"alpha">(d)), get<"ax">(d), get<"keys">(d)); }
        timer.step("mpi::synchronizeHalos");

        release(d, "divv", "curlv");
        acquire(d, "ay", "az");
        computeMomentumEnergy<avClean>(Base::groups_.view(), nullptr, d, domain.box());
        timer.step("MomentumAndEnergy");
        pmReader.step();

        if (d.g != 0.0)
        {
            auto groups = mHolder_.computeSpatialGroups(d, domain);
            mHolder_.upsweep(d, domain);
            timer.step("Upsweep");
            pmReader.step();
            mHolder_.traverse(groups, d, domain);
            timer.step("Gravity");
            pmReader.step();
        }

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
