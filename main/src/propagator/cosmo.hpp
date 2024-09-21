/*
 * MIT License
 *
 * Copyright (c) 2022 CSCS, ETH Zurich
 *               2022 University of Basel
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
 * @brief A Propagator class for plain N-body, computing only gravitational interactions
 *
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

#pragma once

#include <variant>

#include "cstone/fields/field_get.hpp"
#include "sph/particles_data.hpp"
#include "sph/positions.hpp"
// #include "sph/timestep.hpp"

#include "ipropagator.hpp"
#include "gravity_wrapper.hpp"

namespace sphexa
{

using namespace sph;
using util::FieldList;

template<class DomainType, class DataType>
class CosmoProp final : public Propagator<DomainType, DataType>
{
    using Base = Propagator<DomainType, DataType>;
    using Base::timer;

    using KeyType       = typename DataType::KeyType;
    using T             = typename DataType::RealType;
    using HydroType     = typename DataType::HydroData::HydroType;
    using Tmass         = typename DataType::HydroData::Tmass;
    using MultipoleType = ryoanji::CartesianQuadrupole<Tmass>;

    using Acc       = typename DataType::AcceleratorType;
    using MHolder_t = typename cstone::AccelSwitchType<Acc, MultipoleHolderCpu, MultipoleHolderGpu>::template type<
        MultipoleType, DomainType, typename DataType::HydroData>;

    MHolder_t      mHolder_;
    GroupData<Acc> groups_;

    /*! @brief the list of conserved particles fields with values preserved between iterations
     *
     * x, y, z, h and m are automatically considered conserved and must not be specified in this list
     */
    using ConservedFields = FieldList<"vx", "vy", "vz", "vhx", "vhy", "vhz", "du_m1", "u">;

    //! @brief the list of dependent particle fields, these may be used as scratch space during domain sync
    using DependentFields =
        FieldList<"rho", "c", "p", "ax", "ay", "az", "du", "c11", "c12", "c13", "c22", "c23", "c33", "nc">;
    bool halfStepKickNeeded;

public:
    CosmoProp(std::ostream& output, size_t rank)
        : Base(output, rank)
        , halfStepKickNeeded(true)
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

        //! @brief Fields accessed in domain sync are not part of extensible lists.
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
        domain.syncGrav(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                        get<ConservedFields>(d), get<DependentFields>(d));
    }

    void computeGravity(DomainType& domain, DataType& simData)
    {
        auto& d = simData.hydro;

        mHolder_.upsweep(d, domain);
        timer.step("Upsweep");
        mHolder_.traverse(d, domain);
        timer.step("Gravity");

        auto stats = mHolder_.readStats();

        if (domain.startIndex() == 0 && cstone::HaveGpu<typename DataType::AcceleratorType>{})
        {
            //            size_t n = last - first;
            //            std::cout << "numP2P " << stats[0] / n << " maxP2P " << stats[1] << " numM2P " << stats[2] / n
            //            << " maxM2P "
            //                      << stats[3] << std::endl;
        }
    }
    //
    //    void computeAccelerations(size_t first, size_t last, Dataset& d)
    //    {
    //        fill(get<"ax">(d), first, last, HydroType(0));
    //        fill(get<"ay">(d), first, last, HydroType(0));
    //        fill(get<"az">(d), first, last, HydroType(0));
    //
    //        computeGravity(domain, simData);
    //
    //        /* zero other hydro accelerations and compute those */
    //        /* ... */
    //    }
    //
    //    void GravKick(DomainType& domain, double t, double dt)
    //    {
    //        auto dt_cosmo = dt;
    //
    // #pragma omp parallel for schedule(static)
    //        for (size_t i = startIndex; i < endIndex; i++)
    //        {
    //            cstone::Vec3<T> V{d.vx[i], d.vy[i], d.vz[i]};
    //            cstone::Vec3<T> A{d.ax[i], d.ay[i], d.az[i]};
    //
    //            V += A * dt_cosmo;
    //
    //            util::tie(d.vx[i], d.vy[i], d.vz[i]) = util::tie(V[0], V[1], V[2]);
    //        }
    //    }

    //    void GravDrift(DomainType& domain, double t, double dt)
    //    {
    //        auto dt_cosmo = dt;
    //
    // #pragma omp parallel for schedule(static)
    //        for (size_t i = startIndex; i < endIndex; i++)
    //        {
    //            cstone::Vec3<T> X{d.x[i], d.y[i], d.z[i]};
    //            cstone::Vec3<T> V{d.vx[i], d.vy[i], d.vz[i]};
    //
    //            util::tie(d.x_m1[i], d.y_m1[i], d.z_m1[i]) = util::tie(d.x[i], d.y[i], d.z[i]);
    //
    //            X += V * dt_cosmo;
    //            X = cstone::putInBox(X, domain.box()); /* Maybe we do this after all drifts? */
    //
    //            util::tie(d.x[i], d.y[i], d.z[i]) = util::tie(X[0], X[1], X[2]);
    //        }
    //    }
    //

    void gravity_kick(DomainType& domain, DataType& simData, const double dt)
    {
        // Drift is the same for Gravity and SPH
        // Instead of vx use x_m1 * mindt
        auto& d = simData.hydro;

        const size_t first           = domain.startIndex();
        const size_t last            = domain.endIndex();
        const auto   dt_cosmological = dt;

#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; i++)
        {
            cstone::Vec3<T> V{d.vx[i], d.vy[i], d.vz[i]};
            //            util::tie(d.vx_m1[i], d.vy_m1[i], d.vz_m1[i]) = util::tie(d.vx[i], d.vy[i], d.vz[i]);

            const cstone::Vec3<T> A{d.agx[i], d.agy[i], d.agz[i]};

            V += A * dt_cosmological;
            util::tie(d.vx[i], d.vy[i], d.vz[i]) = util::tie(V[0], V[1], V[2]);
        }
    }

    void hydro_force_kick(DomainType& domain, DataType& simData, const double dt)
    {
        // Drift is the same for Gravity and SPH
        auto& d = simData.hydro;

        const size_t first           = domain.startIndex();
        const size_t last            = domain.endIndex();
        const auto   dt_cosmological = dt;

#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; i++)
        {
            cstone::Vec3<T>       V{d.vx[i], d.vy[i], d.vz[i]};
            const cstone::Vec3<T> A{d.ax[i], d.ay[i], d.az[i]};

            V += A * dt_cosmological;
            util::tie(d.vx[i], d.vy[i], d.vz[i]) = util::tie(V[0], V[1], V[2]);
        }
    }

    // Adapt for groups
    void kick(DomainType& domain, DataType& simData, const double dt)
    {
        gravity_kick(domain, simData, dt);
        hydro_force_kick(domain, simData.hydro, dt);
    }

    // Adapt for groups
    void drift(DomainType& domain, DataType& simData, const double dt)
    {
        // Drift is the same for Gravity and SPH
        auto& d = simData.hydro;

        const size_t first    = domain.startIndex();
        const size_t last     = domain.endIndex();
        const auto   dt_cosmo = dt;

#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; i++)
        {
            cstone::Vec3<T> X{d.x[i], d.y[i], d.z[i]};
            cstone::Vec3<T> V{d.vx[i], d.vy[i], d.vz[i]};

            util::tie(d.x_m1[i], d.y_m1[i], d.z_m1[i]) = util::tie(d.x[i], d.y[i], d.z[i]);

            X += V * dt_cosmo;
            X = cstone::putInBox(X, domain.box()); /* Maybe we do this after all drifts? */

            util::tie(d.x[i], d.y[i], d.z[i]) = util::tie(X[0], X[1], X[2]);
            // Integrate U using drift
        }
    }

    void step(DomainType& domain, DataType& simData)
    {
        timer.start();
        sync(domain, simData);
        timer.step("domain::sync");

        auto& d = simData.hydro;
        d.resize(domain.nParticlesWithHalos());
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        transferToHost(d, first, first + 1, {"m"});
        fill(get<"m">(d), 0, first, d.m[first]);
        fill(get<"m">(d), last, domain.nParticlesWithHalos(), d.m[first]);

        if (halfStepKickNeeded)
        {
            //            computeAccelerations(first, last, d);
            timer.step("Accelerations");
            halfStepKickNeeded = false;

            d.minDtCourant = INFINITY;
            computeTimestep(first, last, d);
            timer.step("Timestep");
        }

        //        Kick(domain, d.minDt / 2);
        timer.step("Kick");

        //        Drift(domain, d.minDt);
        timer.step("Drift");

        //        computeAccelerations(first, last, d);
        timer.step("Accelerations");

        d.minDtCourant = INFINITY;
        computeTimestep(first, last, d);
        timer.step("Timestep");

        //        Kick(domain, d.minDt / 2);
        timer.step("Kick");

        timer.step("Leapfrog-Integration");

        // computePositions(first, last, d, domain.box());
        // timer.step("UpdateQuantities");
    }

    void computeForces(DomainType& domain, DataType& simData) override
    {
        timer.start();

        sync(domain, simData);
        timer.step("domain::sync");

        auto& d = simData.hydro;
        d.resize(domain.nParticlesWithHalos());
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        transferToHost(d, first, first + 1, {"m"});
        fill(get<"m">(d), 0, first, d.m[first]);
        fill(get<"m">(d), last, domain.nParticlesWithHalos(), d.m[first]);

        resizeNeighbors(d, domain.nParticles() * d.ngmax);
        findNeighborsSfc(first, last, d, domain.box());
        computeGroups(first, last, d, domain.box(), groups_);
        timer.step("FindNeighbors");

        computeDensity(groups_.view(), d, domain.box());
        timer.step("Density");
        computeEOS_HydroStd(first, last, d);
        timer.step("EquationOfState");

        domain.exchangeHalos(get<"vx", "vy", "vz", "rho", "p", "c">(d), get<"ax">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        computeIAD(groups_.view(), d, domain.box());
        timer.step("IAD");

        domain.exchangeHalos(get<"c11", "c12", "c13", "c22", "c23", "c33">(d), get<"ax">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        computeMomentumEnergySTD(groups_.view(), d, domain.box());
        timer.step("MomentumEnergyIAD");
//Release c11-c13; Acquire agx..agz
        //        if (d.g != 0.0)
        //        {
        //            auto groups = mHolder_.computeSpatialGroups(d, domain);
        //            mHolder_.upsweep(d, domain);
        //            timer.step("Upsweep");
        //            //Add the gravity accelerations to another field
        //            mHolder_.traverse(groups, d, domain);
        //            timer.step("Gravity");
        //        }
    }
    void integrateToFullStep(DomainType& domain, DataType& simData) override
    {
        kick(domain, simData, simData.hydro.minDt / 2.);
    }
    void integrate(DomainType& domain, DataType& simData) override
    {
        auto&  d     = simData.hydro;
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        computeTimestep(first, last, d);
        timer.step("Timestep");

        drift(domain, simData, simData.hydro.minDt);
        kick(domain, simData, simData.hydro.minDt / 2.);
        // Speichere die pred. velocities in vx, vy, vz. Den Kick aber nur mit x - x_m1
        predict_velocities(domain, simData, simData.hydro.minDt / 2.);

        updateSmoothingLength(groups_.view(), d);
        timer.step("UpdateQuantities");
    }
};

} // namespace sphexa
