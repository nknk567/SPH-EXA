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
 * @brief Propagator for
 *
 * @author Noah Kubli <noah.kubli@me.com>
 */

#pragma once

#include <variant>

#include "cstone/fields/particles_get.hpp"
#include "sph/particles_data.hpp"
#include "sph/sph.hpp"

#include "cooling/cooler.hpp"
#include "cooling/eos_cooling.hpp"
#include "ipropagator.hpp"
#include "gravity_wrapper.hpp"

namespace sphexa
{

using namespace sph;
using cstone::FieldList;

template<class DomainType, class DataType>
class CloudProp final : public HydroProp<DomainType, DataType>
{
    using Base = HydroProp<DomainType, DataType>;
    using Base::timer;

    using T       = typename DataType::RealType;
    using KeyType = typename DataType::KeyType;

    cooling::Cooler<T> cooling_data;

    /*! @brief the list of conserved particles fields with values preserved between iterations
     *
     * x, y, z, h and m are automatically considered conserved and must not be specified in this list
     */
    using ConservedFields = FieldList<"u", "vx", "vy", "vz", "x_m1", "y_m1", "z_m1", "du_m1", "alpha", "soft">;

    //! @brief the list of dependent particle fields, these may be used as scratch space during domain sync
    using DependentFields =
        FieldList<"rho", "p", "c", "ax", "ay", "az", "du", "c11", "c12", "c13", "c22", "c23", "c33", "nc">;

    using CoolingFields =
        FieldList<"HI_fraction", "HII_fraction", "HM_fraction", "HeI_fraction", "HeII_fraction", "HeIII_fraction",
                  "H2I_fraction", "H2II_fraction", "DI_fraction", "DII_fraction", "HDI_fraction", "e_fraction",
                  "metal_fraction", "volumetric_heating_rate", "specific_heating_rate", "RT_heating_rate",
                  "RT_HI_ionization_rate", "RT_HeI_ionization_rate", "RT_HeII_ionization_rate",
                  "RT_H2_dissociation_rate", "H2_self_shielding_length">;

public:
    CloudProp(std::ostream& output, size_t rank)
        : Base(output, rank)
    {

        constexpr float ms_sim = 1e8; // 1e9;//1e16;
        constexpr float kp_sim = 1.0; // 1.0;//46400.;

        std::map<std::string, std::any> grackleOptions;
        grackleOptions["use_grackle"]            = 1;
        grackleOptions["with_radiative_cooling"] = 1;
        grackleOptions["primordial_chemistry"]   = 3;
        grackleOptions["dust_chemistry"]         = 0;
        grackleOptions["metal_cooling"]          = 0;
        grackleOptions["UVbackground"]           = 1;
        cooling_data.init(ms_sim, kp_sim, 0, grackleOptions, std::nullopt);
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
        std::apply([&simData](auto... f) { simData.chem.setConserved(f.value...); }, make_tuple(CoolingFields{}));

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
            std::cout << "sizes: " << get<"x">(d).size() << "\t" << get<"HI_fraction">(simData.chem).size()
                      << std::endl;
            domain.syncGrav(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                            std::tuple_cat(get<ConservedFields>(d), get<CoolingFields>(simData.chem)),
                            get<DependentFields>(d));
        }
        else
        {
            domain.sync(
                get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d),
                std::tuple_cat(std::tie(get<"m">(d)), get<ConservedFields>(d), get<CoolingFields>(simData.chem)),
                get<DependentFields>(d));
        }
        d.treeView = domain.octreeProperties();
    }

    void computeForces(DomainType& domain, DataType& simData)
    {
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();
        auto&  d     = simData.hydro;

        resizeNeighbors(d, domain.nParticles() * d.ngmax);
        findNeighborsSfc(first, last, d, domain.box());
        timer.step("FindNeighbors");

        computeDensity(first, last, d, domain.box());
        timer.step("Density");
        // computeEOS_HydroStd(first, last, d);
        eos_cooling(first, last, d, simData.chem, cooling_data);
        timer.step("EquationOfState");

        domain.exchangeHalos(get<"vx", "vy", "vz", "rho", "p", "c">(d), get<"ax">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        computeIAD(first, last, d, domain.box());
        timer.step("IAD");

        domain.exchangeHalos(get<"c11", "c12", "c13", "c22", "c23", "c33">(d), get<"ax">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        computeMomentumEnergySTD(first, last, d, domain.box());
        timer.step("MomentumEnergyIAD");

        if (d.g != 0.0)
        {
            Base::mHolder_.upsweep(d, domain);
            timer.step("Upsweep");
            Base::mHolder_.traverse(d, domain);
            timer.step("Gravity");
        }
    }

    void prepareSystem(DomainType& domain, DataType& simData) override
    {
        auto& d = simData.hydro;
        timer.start();

        sync(domain, simData);
        domain.exchangeHalos(std::tie(get<"m">(d)), get<"ax">(d), get<"ay">(d));

        d.resize(domain.nParticlesWithHalos());
        computeForces(domain, simData);
    }

    void relaxSystem(DomainType& domain, DataType& simData) override
    {
        size_t step = 0;
        while (1)
        {
            auto& d = simData.hydro;
            timer.start();

            sync(domain, simData);
            // halo exchange for masses, allows for particles with variable masses
            domain.exchangeHalos(std::tie(get<"m">(d)), get<"ax">(d), get<"ay">(d));
            timer.step("domain::sync");

            d.resize(domain.nParticlesWithHalos());
            std::cout << get<"u">(d)[0] << std::endl;
            computeForces(domain, simData);

            size_t first = domain.startIndex();
            size_t last  = domain.endIndex();

            computeTimestep(first, last, d);
            timer.step("Timestep");
            // Friction
            const T friction_time{10. * d.minDt};
            T       max_fric{0.};
            T       fric_tot{0.};
            for (size_t i = first; i < last; i++)
            {
                T fric_x = -get<"vx">(d)[i] / friction_time;
                T fric_y = -get<"vy">(d)[i] / friction_time;
                T fric_z = -get<"vz">(d)[i] / friction_time;
                T m_fric{std::max({fric_x, fric_y, fric_z})};
                if (m_fric > max_fric) max_fric = m_fric;
                get<"ax">(d)[i] += fric_x;
                get<"ay">(d)[i] += fric_y;
                get<"az">(d)[i] += fric_z;
                fric_tot += std::abs(fric_x) + std::abs(fric_y) + std::abs(fric_z);
            }
            fric_tot /= d.numParticlesGlobal;
            if (step > 100 && fric_tot < 0.018) break;
            std::cout << "fric_tot " << fric_tot << std::endl;
            computePositions(first, last, d, domain.box());
            timer.step("UpdateQuantities");
            updateSmoothingLength(first, last, d);
            timer.step("UpdateSmoothingLength");

            timer.stop();
            step++;
        }
    }

    void step(DomainType& domain, DataType& simData) override
    {
        auto& d = simData.hydro;
        timer.start();

        sync(domain, simData);
        // halo exchange for masses, allows for particles with variable masses
        domain.exchangeHalos(std::tie(get<"m">(d)), get<"ax">(d), get<"ay">(d));
        timer.step("domain::sync");

        d.resize(domain.nParticlesWithHalos());
        std::cout << get<"u">(d)[0] << std::endl;
        fill(get<"soft">(d), 0, domain.nParticlesWithHalos(), 0.05);
        computeForces(domain, simData);

        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        computeTimestep_cool(first, last, d, cooling_data, simData.chem);
        timer.step("Timestep");

#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; i++)
        {
            //bool haveMui = !d.mui.empty();
            //T    cv      = idealGasCv(haveMui ? d.mui[i] : d.muiConst, d.gamma);

            // T u_old  = cv * d.temp[i];
            // T u_cool = u_old;
            T u_old  = d.u[i];
            T u_cool = d.u[i];
            cooling_data.cool_particle(
                d.minDt, d.rho[i], u_cool, get<"HI_fraction">(simData.chem)[i], get<"HII_fraction">(simData.chem)[i],
                get<"HM_fraction">(simData.chem)[i], get<"HeI_fraction">(simData.chem)[i],
                get<"HeII_fraction">(simData.chem)[i], get<"HeIII_fraction">(simData.chem)[i],
                get<"H2I_fraction">(simData.chem)[i], get<"H2II_fraction">(simData.chem)[i],
                get<"DI_fraction">(simData.chem)[i], get<"DII_fraction">(simData.chem)[i],
                get<"HDI_fraction">(simData.chem)[i], get<"e_fraction">(simData.chem)[i],
                get<"metal_fraction">(simData.chem)[i], get<"volumetric_heating_rate">(simData.chem)[i],
                get<"specific_heating_rate">(simData.chem)[i], get<"RT_heating_rate">(simData.chem)[i],
                get<"RT_HI_ionization_rate">(simData.chem)[i], get<"RT_HeI_ionization_rate">(simData.chem)[i],
                get<"RT_HeII_ionization_rate">(simData.chem)[i], get<"RT_H2_dissociation_rate">(simData.chem)[i],
                get<"H2_self_shielding_length">(simData.chem)[i]);
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
        Base::saveFields(writer, first, last, simData, box);
        //should be customized and namespace
        auto&            chem             = simData.chem;
        auto             fieldPointers = chem.data();
        //std::vector<int> outputFields  = chem.outputFieldIndices;

        auto output = [&]()
        {
            for (int i = int(fieldPointers.size()) - 1; i >= 0; --i)
            {
                //int fidx = outputFields[i];
                //if (d.isAllocated(fidx))
                //{
                    //int column = std::find(d.outputFieldIndices.begin(), d.outputFieldIndices.end(), fidx) -
                     //            d.outputFieldIndices.begin();
                    transferToHost(chem, first, last, {chem.fieldNames[i]});
                    std::visit([writer, i, key = chem.fieldNames[i]](auto field)
                               { writer->writeField(key, field->data(), i); },
                               fieldPointers[i]);
                    //outputFields.erase(outputFields.begin() + i);
                //}
            }
        };

        output();
    }


};

} // namespace sphexa
