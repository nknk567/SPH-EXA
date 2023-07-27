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
// #include "sph/sph.hpp"
//  #include "sph/eos.hpp"
#include "sph/particles_data.hpp"
#include "isim_init.hpp"
#include "early_sync.hpp"
#include "grid.hpp"
#include "cooling/init_chemistry.h"

namespace sphexa
{

std::map<std::string, double> cloudConstants()
{
    return {{"gravConstant", 1.},  {"r", 3.},       {"mTotal", 1.},           {"gamma", 5. / 3.},
            {"u0", 0.05},          {"minDt", 1e-4}, {"minDt_m1", 1e-4},       {"mui", 10},
            {"ng0", 100},          {"ngmax", 110},  {"metal_fraction", 1e-6}, {"hydrogen_fraction", 0.76},
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

    const T u_guess{0.2};
    std::fill(d.u.begin(), d.u.end(), u_guess);

    cooling::initChemistryData(chem, d.x.size());

    T totalVolume = 4 * M_PI / 3 * std::pow(constants.at("r"), 3);
    // before the contraction with sqrt(r), the sphere has a constant particle concentration of Ntot / Vtot
    // after shifting particles towards the center by factor sqrt(r), the local concentration becomes
    // c(r) = 2/3 * 1/r * Ntot / Vtot
    T c0 = 2. / 3. * d.numParticlesGlobal / totalVolume;
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

template<typename FT, typename T>
T bisect_monotone(const FT& func, const T y)
{
    T       x       = 0.;
    T       delta   = 2.;
    const T epsilon = 1e-5;
    while (true)
    {
        // std::cout << "bisect: " << std::endl;
        T       x_new   = x + delta;
        const T y_x_new = func(x_new);
        if (std::abs(y - y_x_new) < epsilon)
        {
            x = x_new;
            break;
        }
        if (y_x_new > y)
            delta *= 0.5;
        else
            x = x_new;
    }
    return x;
}

template<class Vector>
void contractRhoProfileCloud(Vector& x, Vector& y, Vector& z)
{
    // truncated 1/r
    auto f = [](double r) {
        if (r <= 1.) return r;
        else return std::sqrt(1./3.) * std::sqrt(2.*r*r*r + 1.);
    };
    // Exponential profile
    /*auto f_1 = [](double y) {
        //return 2. - std::exp(-y) * (y*y + 2.*y + 2.);
        //return std::pow((2. - std::exp(-y) * (y*y + 2.*y + 2.)) * 3., 1./3.);
        const double b = 1.75;
        return std::pow((2./(b*b*b) - std::exp(-b * y) * (b*b*y*y + 2.*b*y + 2.) / (b*b*b)) * 3., 1./3.);
    };*/
    // double x = bisect_monotone(f, 1.0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < x.size(); i++)
    {
        auto radius0 = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

        // multiply coordinates by sqrt(r) to generate a density profile ~ 1/r
        // auto contraction = std::sqrt(radius0);
        //auto new_r       = bisect_monotone(f_1, radius0);
        auto new_r       = f(radius0);
        auto contraction = new_r / radius0;
        x[i] *= contraction;
        y[i] *= contraction;
        z[i] *= contraction;
    }
}

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

    bool initDependent(Dataset& simData) const override
    {
        using T    = typename Dataset::RealType;
        auto& d    = simData.hydro;
        auto& chem = simData.chem;

        // const T        rho_const{3. / (4. * M_PI)};
        std::vector<T> pressure_eq(d.x.size(), 0.);
        for (size_t i = 0; i < d.x.size(); i++)
        {
            const T radius = std::sqrt(d.x[i] * d.x[i] + d.y[i] * d.y[i] + d.z[i] * d.z[i]);
            // pressure_eq[i] = 3. / (8. * M_PI) * (1. - radius * radius);

            //exponential profile
            /*auto p = [](double radius)
            {
                const T b = 1.75;
                // const T e = std::exp(b * radius);
                // const T e2 = std::exp(2. * b * radius);
                const T em   = std::exp(-b * radius);
                const T em2  = std::exp(-2. * b * radius);
                const T eim2 = std::expint(-2 * b * radius);
                const T eim  = std::expint(-b * radius);

                const T t1  = 4. * b * radius * eim2;
                const T t2  = -4. * b * radius * eim;
                const T t3  = em2 * (b * radius + 4.) - 4. * em;
                const T n   = 2. * M_PI * (t1 + t2 + t3);
                const T res = (n / (radius * b * b * b) * (9. / (16. * M_PI * M_PI)));
                return res;
            };*/
           // pressure_eq[i] = p(2.9892) - p(radius);
            //truncated 1/r
            const T r0 = settings_.at("r");
            const double rc = std::sqrt(1. / 3.) * std::sqrt(2. * r0*r0*r0 + 2.);
            const double k = 3. / (4. * M_PI * r0*r0*r0);
            const double B = 2. * M_PI *  k*k * (std::log(rc) - std::log(1.) + 1. / (6. * rc * rc) - 1. / 6.);
            auto p = [&](double radius)
            {
             if (radius <= 1.) return B + 2. * M_PI * k*k / 3. * (1. - radius * radius);
             else return 2. * M_PI *  k*k * (std::log(rc) - std::log(radius) + 1. / (6. * rc * rc) - 1. / (6. * radius * radius));
            };

            pressure_eq[i] = p(radius);

            /*const T ei = std::expint(-2. * b * radius);
            const T e = std::exp(-2 * b * radius);
            const T first_term = -2 * b * ei;
            const T second_term = -e * (b * radius + 4.0) / (2. * radius);
            const T pressure = 0. - (first_term + second_term) / (b*b*b);
            pressure_eq[i] = pressure * 3. / (4 * M_PI);*/

            /*const T r_eps = std::max(radius, 1e-5);
            pressure_eq[i] = 1. / (2. * M_PI) * std::log(1. / r_eps);*/
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
        bool good = true;
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
            if (std::abs(relation - 1.) > 1e-5) good = false;
            std::cout << "u: " << d.u[i] << std::endl;
            std::cout << "relation " << relation << std::endl;
            // std::cout << d.rho[i] << std::endl;
        }
        std::cout << "status " << good << std::endl;
        return good;
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

        contractRhoProfileCloud(d.x, d.y, d.z);
        syncCoords<KeyType>(rank, numRanks, numParticlesGlobal, d.x, d.y, d.z, globalBox);

        d.resize(d.x.size());

        settings_["numParticlesGlobal"] = double(numParticlesGlobal);
        std::cout << "numParticlesGlobal: " << numParticlesGlobal << std::endl;
        BuiltinWriter attributeSetter(settings_);
        d.loadOrStoreAttributes(&attributeSetter);

        initCloudFields(d, simData.chem, settings_, globalBox);

        return globalBox;
    }

    const std::map<std::string, double>& constants() const override { return settings_; }
};

} // namespace sphexa
