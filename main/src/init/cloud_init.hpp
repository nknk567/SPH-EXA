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
 * @brief Evrard collapse initialization
 *
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

#pragma once

#include <map>

#include "cstone/sfc/box.hpp"
#include "cstone/tree/continuum.hpp"
#include "sph/sph.hpp"
// #include "sph/eos.hpp"
#include "sph/particles_data.hpp"
#include "isim_init.hpp"
#include "early_sync.hpp"
#include "grid.hpp"
#include "cooling/init_chemistry.h"

namespace sphexa
{

std::map<std::string, double> cloudConstants()
{
    return {{"gravConstant", 1.},  {"r", 1.},       {"mTotal", 1.},           {"gamma", 5. / 3.},
            {"u0", 0.05},          {"minDt", 1e-4}, {"minDt_m1", 1e-4},       {"mui", 10},
            {"ng0", 100},          {"ngmax", 150},  {"metal_fraction", 1e-6}, {"hydrogen_fraction", 0.76},
            {"d_to_h_ratio", 1e-5}};
}

template<class Dataset, typename ChemData>
void initCloudFields(Dataset& d, ChemData& chem, const std::map<std::string, double>& constants,
                     const cstone::Box<typename Dataset::RealType>& box)
{
    using T = typename Dataset::RealType;

    double mPart = constants.at("mTotal") / d.numParticlesGlobal;

    std::fill(d.m.begin(), d.m.end(), mPart);
    std::fill(d.du_m1.begin(), d.du_m1.end(), 0.0);
    std::fill(d.mui.begin(), d.mui.end(), d.muiConst);
    std::fill(d.alpha.begin(), d.alpha.end(), d.alphamin);

    std::fill(d.vx.begin(), d.vx.end(), 0.0);
    std::fill(d.vy.begin(), d.vy.end(), 0.0);
    std::fill(d.vz.begin(), d.vz.end(), 0.0);

    std::fill(d.x_m1.begin(), d.x_m1.end(), 0.0);
    std::fill(d.y_m1.begin(), d.y_m1.end(), 0.0);
    std::fill(d.z_m1.begin(), d.z_m1.end(), 0.0);

    const T u_guess{2.87};
    std::fill(d.u.begin(), d.u.end(), u_guess);
    /*auto cv    = sph::idealGasCv(d.muiConst, d.gamma);
    auto temp0 = constants.at("u0") / cv;
    std::fill(d.temp.begin(), d.temp.end(), temp0);*/

    cooling::initChemistryData(chem, d.x.size());
    // Calculate u
    /*const T        rho_const{3. / (4. * M_PI)};
    std::vector<T> pressure_eq(d.x.size(), 0.);
    for (size_t i = 0; i < d.x.size(); i++)
    {
        const T radius = std::sqrt(d.x[i] * d.x[i] + d.y[i] * d.y[i] + d.z[i] * d.z[i]);
        pressure_eq[i] = 3. / (8. * M_PI) * (1. - radius * radius);
    }



    cooling::Cooler<T>              cooling_data;
    constexpr float                 ms_sim = 1e6;
    constexpr float                 kp_sim = 1.0;
    std::map<std::string, std::any> grackleOptions;
    grackleOptions["use_grackle"]            = 1;
    grackleOptions["with_radiative_cooling"] = 0;
    grackleOptions["primordial_chemistry"]   = 3;
    grackleOptions["dust_chemistry"]         = 0;
    grackleOptions["metal_cooling"]          = 1;
    grackleOptions["UVbackground"]           = 1;
    cooling_data.init(ms_sim, kp_sim, 0, grackleOptions, std::nullopt);

    T nden = get<"metal_fraction">(chem)[0] * rho_const / 16.;
    nden += (get<"HI_fraction">(chem)[0] + get<"HII_fraction">(chem)[0] + get<"e_fraction">(chem)[0] +
             (get<"HeI_fraction">(chem)[0] + get<"HeII_fraction">(chem)[0] + get<"HeIII_fraction">(chem)[0]) / 4.) *
            rho_const;
    nden +=
        (get<"HM_fraction">(chem)[0] + (get<"H2I_fraction">(chem)[0] + get<"H2II_fraction">(chem)[0]) / 2.) * rho_const;
    const T mu = rho_const / nden;
    const T u_guess{2.87};
    //Calculate density
    //resizeNeighbors(d, d.x.size() * d.ngmax);
    //sph::findNeighborsSfc(size_t(0), d.x.size(), d, box);
    //sph::computeDensity<T, Dataset>(size_t(0), d.x.size(), d, box);
    for (size_t i = 0; i < d.x.size(); i++)
    {
        T rho{rho_const};
        T u_cool{u_guess};
        const T pressure = cooling_data.pressure(
            rho, u_cool, get<"HI_fraction">(chem)[i], get<"HII_fraction">(chem)[i],
            get<"HM_fraction">(chem)[i], get<"HeI_fraction">(chem)[i],
            get<"HeII_fraction">(chem)[i], get<"HeIII_fraction">(chem)[i],
            get<"H2I_fraction">(chem)[i], get<"H2II_fraction">(chem)[i],
            get<"DI_fraction">(chem)[i], get<"DII_fraction">(chem)[i],
            get<"HDI_fraction">(chem)[i], get<"e_fraction">(chem)[i],
            get<"metal_fraction">(chem)[i], get<"volumetric_heating_rate">(chem)[i],
            get<"specific_heating_rate">(chem)[i], get<"RT_heating_rate">(chem)[i],
            get<"RT_HI_ionization_rate">(chem)[i], get<"RT_HeI_ionization_rate">(chem)[i],
            get<"RT_HeII_ionization_rate">(chem)[i], get<"RT_H2_dissociation_rate">(chem)[i],
            get<"H2_self_shielding_length">(chem)[i]);
        const T relation{pressure_eq[i] / pressure};
        d.u[i] = u_guess * relation;
        //std::cout << d.u[i] << std::endl;
        //std::cout << d.rho[i] << std::endl;
    }*/

    T totalVolume = 4 * M_PI / 3 * std::pow(constants.at("r"), 3);
    // before the contraction with sqrt(r), the sphere has a constant particle concentration of Ntot / Vtot
    // after shifting particles towards the center by factor sqrt(r), the local concentration becomes
    // c(r) = 2/3 * 1/r * Ntot / Vtot
    T c0 = /*2. / 3. * */ d.numParticlesGlobal / totalVolume;
    std::cout << d.numParticlesGlobal << "numParticlesGlobal" << std::endl;
    std::cout << c0 << "c0" << std::endl;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < d.x.size(); i++)
    {
        T radius        = std::sqrt((d.x[i] * d.x[i]) + (d.y[i] * d.y[i]) + (d.z[i] * d.z[i]));
        T concentration = c0 / radius;
        d.h[i]          = std::cbrt(3 / (4 * M_PI) * d.ng0 / concentration) * 0.5;
        std::cout << "h: " << d.h[i] << std::endl; // is it infinity?
    }
}

/*template<class Vector>
void contractRhoProfile(Vector& x, Vector& y, Vector& z)
{
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < x.size(); i++)
    {
        auto radius0 = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

        // multiply coordinates by sqrt(r) to generate a density profile ~ 1/r
        auto contraction = std::sqrt(radius0);
        x[i] *= contraction;
        y[i] *= contraction;
        z[i] *= contraction;
    }
}*/

/*//! @brief Estimate SFC partition of the Evrard sphere based on approximate continuum particle counts
template<class KeyType, class T>
std::tuple<KeyType, KeyType> estimateEvrardSfcPartition(size_t cbrtNumPart, const cstone::Box<T>& box, int rank,
                                                        int numRanks)
{
    size_t numParticlesGlobal = 0.523 * cbrtNumPart * cbrtNumPart * cbrtNumPart;
    T      r                  = box.xmax();

    double   eps        = 2.0 * r / (1u << cstone::maxTreeLevel<KeyType>{});
    unsigned bucketSize = numParticlesGlobal / (100 * numRanks);

    auto oneOverR = [numParticlesGlobal, r, eps](T x, T y, T z)
    {
        T radius = std::max(std::sqrt(norm2(cstone::Vec3<T>{x, y, z})), eps);
        if (radius > r) { return 0.0; }
        else { return T(numParticlesGlobal) / (2 * M_PI * radius); }
    };

    auto [tree, counts]            = cstone::computeContinuumCsarray<KeyType>(oneOverR, box, bucketSize);
    cstone::SpaceCurveAssignment a = cstone::singleRangeSfcSplit(counts, numRanks);

    return {tree[a.firstNodeIdx(rank)], tree[a.lastNodeIdx(rank)]};
}*/

template<class Dataset>
class CloudGlassSphere : public ISimInitializer<Dataset>
{
    std::string glassBlock;
    using Base = ISimInitializer<Dataset>;
    using Base::settings_;

public:
    explicit CloudGlassSphere(std::string initBlock)
        : glassBlock(std::move(initBlock))
    {
        Base::updateSettings(cloudConstants());
    }

    void initDependent(Dataset& simData) const override
    {
        using T    = typename Dataset::RealType;
        auto& d    = simData.hydro;
        auto& chem = simData.chem;

        const T        rho_const{3. / (4. * M_PI)};
        std::vector<T> pressure_eq(d.x.size(), 0.);
        for (size_t i = 0; i < d.x.size(); i++)
        {
            const T radius = std::sqrt(d.x[i] * d.x[i] + d.y[i] * d.y[i] + d.z[i] * d.z[i]);
            pressure_eq[i] = 3. / (8. * M_PI) * (1. - radius * radius);
        }

        /*cooling::Cooler<T>              cooling_data;
        constexpr float                 ms_sim = 1e6;
        constexpr float                 kp_sim = 1.0;
        std::map<std::string, std::any> grackleOptions;
        grackleOptions["use_grackle"]            = 1;
        grackleOptions["with_radiative_cooling"] = 0;
        grackleOptions["primordial_chemistry"]   = 3;
        grackleOptions["dust_chemistry"]         = 0;
        grackleOptions["metal_cooling"]          = 1;
        grackleOptions["UVbackground"]           = 1;
        cooling_data.init(ms_sim, kp_sim, 0, grackleOptions, std::nullopt);*/

        /*T nden = get<"metal_fraction">(chem)[0] * rho_const / 16.;
        nden += (get<"HI_fraction">(chem)[0] + get<"HII_fraction">(chem)[0] + get<"e_fraction">(chem)[0] +
                 (get<"HeI_fraction">(chem)[0] + get<"HeII_fraction">(chem)[0] + get<"HeIII_fraction">(chem)[0]) / 4.) *
                rho_const;
        nden += (get<"HM_fraction">(chem)[0] + (get<"H2I_fraction">(chem)[0] + get<"H2II_fraction">(chem)[0]) / 2.) *
                rho_const;
        const T mu = rho_const / nden;*/
        // const T u_guess{2.87};
        //  Calculate density
        //  resizeNeighbors(d, d.x.size() * d.ngmax);
        //  sph::findNeighborsSfc(size_t(0), d.x.size(), d, box);
        //  sph::computeDensity<T, Dataset>(size_t(0), d.x.size(), d, box);
        for (size_t i = 0; i < d.x.size(); i++)
        {
            T rho{d.rho[i]};
            T u_cool{d.u[i]};
            T pressure{d.p[i]};
            /*const T pressure = cooling_data.pressure(
                rho, u_cool, get<"HI_fraction">(chem)[i], get<"HII_fraction">(chem)[i], get<"HM_fraction">(chem)[i],
                get<"HeI_fraction">(chem)[i], get<"HeII_fraction">(chem)[i], get<"HeIII_fraction">(chem)[i],
                get<"H2I_fraction">(chem)[i], get<"H2II_fraction">(chem)[i], get<"DI_fraction">(chem)[i],
                get<"DII_fraction">(chem)[i], get<"HDI_fraction">(chem)[i], get<"e_fraction">(chem)[i],
                get<"metal_fraction">(chem)[i], get<"volumetric_heating_rate">(chem)[i],
                get<"specific_heating_rate">(chem)[i], get<"RT_heating_rate">(chem)[i],
                get<"RT_HI_ionization_rate">(chem)[i], get<"RT_HeI_ionization_rate">(chem)[i],
                get<"RT_HeII_ionization_rate">(chem)[i], get<"RT_H2_dissociation_rate">(chem)[i],
                get<"H2_self_shielding_length">(chem)[i]);*/
            const T relation{pressure_eq[i] / pressure};
            d.u[i] = u_cool * relation;
            std::cout << d.u[i] << std::endl;
            std::cout << d.rho[i] << std::endl;
        }
    }
    cstone::Box<typename Dataset::RealType> init(int rank, int numRanks, size_t cbrtNumPart,
                                                 Dataset& simData) const override
    {
        auto& d       = simData.hydro;
        using KeyType = typename Dataset::KeyType;
        using T       = typename Dataset::RealType;

        std::vector<T> xBlock, yBlock, zBlock;
        fileutils::readTemplateBlock(glassBlock, xBlock, yBlock, zBlock);
        size_t blockSize = xBlock.size();

        int               multi1D      = std::rint(cbrtNumPart / std::cbrt(blockSize));
        cstone::Vec3<int> multiplicity = {multi1D, multi1D, multi1D};

        T              r = settings_.at("r");
        cstone::Box<T> globalBox(-r, r, cstone::BoundaryType::open);

        auto [keyStart, keyEnd] = equiDistantSfcSegments<KeyType>(rank, numRanks, 100);
        assembleCuboid<T>(keyStart, keyEnd, globalBox, multiplicity, xBlock, yBlock, zBlock, d.x, d.y, d.z);
        cutSphere(r, d.x, d.y, d.z);

        size_t numParticlesGlobal = d.x.size();
        MPI_Allreduce(MPI_IN_PLACE, &numParticlesGlobal, 1, MpiType<size_t>{}, MPI_SUM, simData.comm);

        // contractRhoProfile(d.x, d.y, d.z);
        syncCoords<KeyType>(rank, numRanks, numParticlesGlobal, d.x, d.y, d.z, globalBox);

        d.resize(d.x.size());

        settings_["numParticlesGlobal"] = double(numParticlesGlobal);
        BuiltinWriter attributeSetter(settings_);
        d.loadOrStoreAttributes(&attributeSetter);

        initCloudFields(d, simData.chem, settings_, globalBox);

        return globalBox;
    }

    const std::map<std::string, double>& constants() const override { return settings_; }
};

} // namespace sphexa
