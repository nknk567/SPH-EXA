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
#include "cosmo/cosmology_data.hpp"
#include "cosmo/factory.hpp"
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

    // Now just load a default universe
    std::unique_ptr<cosmo::Cosmology<double>> cosmology_ptr;

    /*! @brief the list of conserved particles fields with values preserved between iterations
     *
     * x, y, z, h and m are automatically considered conserved and must not be specified in this list
     */
    using ConservedFields = FieldList<"vx", "vy", "vz", "vhx", "vhy", "vhz", "u">;

    //! @brief the list of dependent particle fields, these may be used as scratch space during domain sync
    using DependentFields =
        FieldList<"rho", "c", "p", "ax", "ay", "az", "agx", "agy", "agz", "c11", "c12", "c13", "du", "nc">;
    bool halfStepKickNeeded;

public:
    CosmoProp(std::ostream& output, size_t rank)
        : Base(output, rank)
        , halfStepKickNeeded(true)
        , cosmology_ptr{cosmo::cosmologyFactory(67.66, 0.3111, 0., 0.6889)}
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
    // Switch in this function-> CPU / GPU
    template<typename T_field_t, typename T_d_field_t>
    void euler_forward(T_field_t& x, const T_d_field_t& x_dot, const size_t first, const size_t last,
                       const size_t delta_t)
    {
#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; i++)
        {
            x[i] += x_dot[i] * delta_t;
        }
    }

    template<typename T_field_t, typename T_d_field_t, typename T>
    void euler_forward(T_field_t& x, T_field_t& y, T_field_t& z, T_d_field_t& x_dot, T_d_field_t& y_dot,
                       T_d_field_t& z_dot, const size_t first, const size_t last, const T delta_t,
                       const cstone::Box<T>& box)
    {
#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; i++)
        {
            x[i] += x_dot[i] * delta_t;
            y[i] += y_dot[i] * delta_t;
            z[i] += z_dot[i] * delta_t;

            const auto X                = cstone::putInBox(cstone::Vec3<T>{x[i], y[i], z[i]}, box);
            util::tie(x[i], y[i], z[i]) = util::tie(X[0], X[1], X[2]);
        }
    }

    void gravity_kick(DomainType& domain, DataType& simData, const double dt)
    {
        const size_t first           = domain.startIndex();
        const size_t last            = domain.endIndex();
        auto&        d               = simData.hydro;
        const auto   dt_cosmological = cosmology_ptr->kickTimeCorrection(d.ttot, dt);

        euler_forward(get<"vhx">(d), get<"agx">(d), first, last, dt_cosmological);
        euler_forward(get<"vhy">(d), get<"agy">(d), first, last, dt_cosmological);
        euler_forward(get<"vhz">(d), get<"agz">(d), first, last, dt_cosmological);
    }

    void hydro_force_kick(DomainType& domain, DataType& simData, const double dt)
    {
        const size_t first           = domain.startIndex();
        const size_t last            = domain.endIndex();
        auto&        d               = simData.hydro;
        const auto   dt_cosmological = cosmology_ptr->kickTimeCorrectionSPH(d.ttot, dt, 5. / 3.);

        euler_forward(get<"vhx">(d), get<"ax">(d), first, last, dt_cosmological);
        euler_forward(get<"vhy">(d), get<"ay">(d), first, last, dt_cosmological);
        euler_forward(get<"vhz">(d), get<"az">(d), first, last, dt_cosmological);
    }

    void internal_energy_kick(DomainType& domain, DataType& simData, const double dt)
    {
        const size_t first           = domain.startIndex();
        const size_t last            = domain.endIndex();
        auto&        d               = simData.hydro;
        const auto   dt_cosmological = cosmology_ptr->driftTimeCorrection(d.ttot, dt);

        euler_forward(get<"u">(d), get<"du">(d), first, last, dt_cosmological);
    }

    void predict_velocities(DomainType& domain, DataType& simData, const double dt)
    {
        const size_t first        = domain.startIndex();
        const size_t last         = domain.endIndex();
        auto&        d            = simData.hydro;
        const auto   dt_corr_grav = cosmology_ptr->kickTimeCorrection(d.ttot, dt);
        const auto   dt_corr_sph  = cosmology_ptr->kickTimeCorrectionSPH(d.ttot, dt, 5. / 3.);

        euler_forward(get<"vx">(d), get<"agx">(d), first, last, dt_corr_grav);
        euler_forward(get<"vy">(d), get<"agy">(d), first, last, dt_corr_grav);
        euler_forward(get<"vz">(d), get<"agz">(d), first, last, dt_corr_grav);
        euler_forward(get<"vx">(d), get<"ax">(d), first, last, dt_corr_sph);
        euler_forward(get<"vy">(d), get<"ay">(d), first, last, dt_corr_sph);
        euler_forward(get<"vz">(d), get<"az">(d), first, last, dt_corr_sph);
    }

    // Adapt for groups
    void drift(DomainType& domain, DataType& simData, const double dt)
    {
        // Drift is the same for Gravity and SPH
        auto& d = simData.hydro;

        const size_t first           = domain.startIndex();
        const size_t last            = domain.endIndex();
        const auto   dt_cosmological = cosmology_ptr->driftTimeCorrection(d.ttot, dt);
        euler_forward(get<"x">(d), get<"y">(d), get<"z">(d), get<"vx">(d), get<"vy">(d), get<"vz">(d), first, last,
                      dt_cosmological, domain.box());
        // #pragma omp parallel for schedule(static)
        //         for (size_t i = first; i < last; i++)
        //         {
        //             cstone::Vec3<T> X{d.x[i], d.y[i], d.z[i]};
        //             cstone::Vec3<T> V{d.vx[i], d.vy[i], d.vz[i]};
        //
        //             util::tie(d.x_m1[i], d.y_m1[i], d.z_m1[i]) = util::tie(d.x[i], d.y[i], d.z[i]);
        //
        //             X += V * dt_cosmo;
        //             X = cstone::putInBox(X, domain.box()); /* Maybe we do this after all drifts? */
        //
        //             util::tie(d.x[i], d.y[i], d.z[i]) = util::tie(X[0], X[1], X[2]);
        //             // Integrate U using drift
        //         }
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

        release(d, "agx", "agy", "agz");
        acquire(d, "c22", "c23", "c33");

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

        release(d, "c22", "c23", "c33");
        acquire(d, "agx", "agy", "agz");

        fill(get<"agx">(d), 0., first, last);
        fill(get<"agy">(d), 0., first, last);
        fill(get<"agz">(d), 0., first, last);

        // Release c11-c13; Acquire agx..agz
        if (d.g != 0.0)
        {
            auto groups = mHolder_.computeSpatialGroups(d, domain);
            mHolder_.upsweep(d, domain);
            timer.step("Upsweep");
            //            //Add the gravity accelerations to another field
            mHolder_.traverse(groups, d, domain);
            timer.step("Gravity");
        }
    }

    //! @brief Implements the first two parts of the kick-drift-kick scheme.
    void integrate(DomainType& domain, DataType& simData) override
    {
        if (simData.hydro.leapfrog_synced) { simData.hydro.leapfrog_synced = false; }
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        computeTimestep(first, last, simData.hydro);
        timer.step("Timestep");

        gravity_kick(domain, simData, simData.hydro.minDt / 2.);
        hydro_force_kick(domain, simData, simData.hydro.minDt / 2.);
        internal_energy_kick(domain, simData, simData.hydro.minDt / 2.);
        drift(domain, simData, simData.hydro.minDt);
        predict_velocities(domain, simData, simData.hydro.minDt);
    }

    //! @brief Implements the final kick in the kick-drift-kick scheme. This kick requires
    //! the new accelerations from computeForces().
    //! Outputs should be written after this step as then the quantities are synchronized.
    //! If one starts from an input file which does not contain both types of velocities
    //! one can set them equal (vx = vhx, ...) and omit this step for the first iteration.
    void integrateToFullStep(DomainType& domain, DataType& simData) override
    {
        if (simData.hydro.leapfrog_synced) { return; }
        gravity_kick(domain, simData, simData.hydro.minDt / 2.);
        hydro_force_kick(domain, simData, simData.hydro.minDt / 2.);
    }
};

} // namespace sphexa
