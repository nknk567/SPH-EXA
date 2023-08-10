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
    return {{"gravConstant", 1.},  {"r", 1.5},      {"mTotal", 1.},           {"gamma", 5. / 3.},
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
        // std::cout << "h: " << d.h[i] << std::endl; // is it infinity?
    }
}

template<class Vector>
void contractRhoProfileCloud(Vector& x, Vector& y, Vector& z)
{
    // truncated 1/r^4
    auto f = [](const double r)
    {
        if (r <= 1.)
            return r;
        else
            return 1. / (4. / 3. - r * r * r / 3.);
    };

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < x.size(); i++)
    {
        const auto radius0 = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

        // multiply coordinates by sqrt(r) to generate the density profile
        auto new_r       = f(radius0);
        auto contraction = new_r / radius0;
        x[i] *= contraction;
        y[i] *= contraction;
        z[i] *= contraction;
    }
}

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

        std::vector<T> pressure_eq(d.x.size(), 0.);
        for (size_t i = 0; i < d.x.size(); i++)
        {
            const T radius = std::sqrt(d.x[i] * d.x[i] + d.y[i] * d.y[i] + d.z[i] * d.z[i]);

            // truncated 1/r^4
            const T      r0 = settings_.at("r");
            const double rc = 1. / (4. / 3. - r0 * r0 * r0 / 3.); // nicht grösser als 4^(1/3)
            const double k  = 3. / (4. * M_PI * r0 * r0 * r0);    // Mtot == 1

            auto Bf = [&](double r)
            { return -4. * M_PI * k * k * (-1. / (6. * std::pow(r, 6.)) + 4. / 15. / std::pow(r, 5.)); };

            const double B     = Bf(rc) - Bf(1.);
            auto         Afunc = [&](double r) { return 2. * M_PI * k * k / 3. * (1. - r * r); };
            auto         p     = [&](double radius)
            {
                if (radius <= 1.)
                    return B + Afunc(radius);
                else
                    return Bf(rc) - Bf(radius);
            };

            pressure_eq[i] = p(radius);
        }

        bool good = true;
        for (size_t i = 0; i < d.x.size(); i++)
        {
            T       rho{d.rho[i]};
            T       u_cool{d.u[i]};
            T       pressure{d.p[i]};
            const T relation{pressure_eq[i] / pressure};
            d.u[i] = u_cool * relation;
            if (std::abs(relation - 1.) > 1e-5) good = false;
        }
        std::cout << "status " << good << std::endl;
        return good;
    }

    void initPressure(Dataset& simData, const int rank, const int numRanks,
                      const cstone::Box<typename Dataset::RealType>& globalBox, const size_t numParticlesGlobal) const
    {
        auto& d       = simData.hydro;
        using KeyType = typename Dataset::KeyType;
        using T       = typename Dataset::RealType;
        using cstone::FieldList;
        using ConservedFields = FieldList<"temp", "vx", "vy", "vz", "x_m1", "y_m1", "z_m1", "du_m1", "u">;

        //! @brief the list of dependent particle fields, these may be used as scratch space during domain sync
        using DependentFields = FieldList<"rho", "p", "c", "du", "nc">;

        using CoolingFields =
            FieldList<"HI_fraction", "HII_fraction", "HM_fraction", "HeI_fraction", "HeII_fraction", "HeIII_fraction",
                      "H2I_fraction", "H2II_fraction", "DI_fraction", "DII_fraction", "HDI_fraction", "e_fraction",
                      "metal_fraction", "volumetric_heating_rate", "specific_heating_rate", "RT_heating_rate",
                      "RT_HI_ionization_rate", "RT_HeI_ionization_rate", "RT_HeII_ionization_rate",
                      "RT_H2_dissociation_rate", "H2_self_shielding_length">;
        d.setConserved("x", "y", "z", "h", "m");
        d.setDependent("keys");
        std::apply([&d](auto... f) { d.setConserved(f.value...); }, make_tuple(ConservedFields{}));
        std::apply([&d](auto... f) { d.setDependent(f.value...); }, make_tuple(DependentFields{}));
        std::apply([&simData](auto... f) { simData.chem.setConserved(f.value...); }, make_tuple(CoolingFields{}));

        d.resize(d.x.size());

        settings_["numParticlesGlobal"] = double(numParticlesGlobal);
        std::cout << "numParticlesGlobal: " << numParticlesGlobal << std::endl;
        BuiltinWriter attributeSetter(settings_);
        d.loadOrStoreAttributes(&attributeSetter);

        initCloudFields(d, simData.chem, settings_, globalBox);

        // Prepare hydrost. eq
        uint64_t bucketSizeFocus = 64;

        uint64_t bucketSize = std::max(bucketSizeFocus, d.numParticlesGlobal / (100 * numRanks));
        cstone::Domain<KeyType, T, cstone::CpuTag> domain(rank, numRanks, bucketSize, bucketSizeFocus, 0.5, globalBox);

        if (rank == 0) std::cout << "nLocalParticles " << get<"HI_fraction">(simData.chem).size() << std::endl;

        auto bl = []()
        {
            std::vector<std::string> ret{"x", "y", "z", "h", "m"};
            for_each_tuple([&ret](auto f) { ret.push_back(f.value); }, make_tuple(ConservedFields{}));
            return ret;
        };

       // transferToDevice(d, 0, d.x.size(), bl());

        domain.syncGrav(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                        std::tuple_cat(get<ConservedFields>(d), get<CoolingFields>(simData.chem)),
                        get<DependentFields>(d));
        d.treeView = domain.octreeProperties();

        d.resize(domain.nParticlesWithHalos());

        size_t                          first = domain.startIndex();
        size_t                          last  = domain.endIndex();
        cooling::Cooler<T>              cooling_data;
        constexpr float                 ms_sim = 1e8;
        constexpr float                 kp_sim = 1.0;
        std::map<std::string, std::any> grackleOptions;
        grackleOptions["use_grackle"]            = 1;
        grackleOptions["with_radiative_cooling"] = 0;
        grackleOptions["primordial_chemistry"]   = 1;
        grackleOptions["dust_chemistry"]         = 0;
        grackleOptions["metal_cooling"]          = 0;
        grackleOptions["UVbackground"]           = 0;

        resizeNeighbors(d, domain.nParticles() * d.ngmax);
        sph::findNeighborsSfc(first, last, d, domain.box());
        sph::computeDensity(first, last, d, domain.box()); // halo exchange rho!!
        cooling_data.init(ms_sim, kp_sim, 0, grackleOptions, std::nullopt);

        auto calculatePressure = [&, &chem = simData.chem]()
        {
#pragma omp parallel for schedule(static)
            for (size_t i = first; i < last; ++i)
            {
                T pressure = cooling_data.pressure(
                    d.rho[i], d.u[i], get<"HI_fraction">(chem)[i], get<"HII_fraction">(chem)[i],
                    get<"HM_fraction">(chem)[i], get<"HeI_fraction">(chem)[i], get<"HeII_fraction">(chem)[i],
                    get<"HeIII_fraction">(chem)[i], get<"H2I_fraction">(chem)[i], get<"H2II_fraction">(chem)[i],
                    get<"DI_fraction">(chem)[i], get<"DII_fraction">(chem)[i], get<"HDI_fraction">(chem)[i],
                    get<"e_fraction">(chem)[i], get<"metal_fraction">(chem)[i], get<"volumetric_heating_rate">(chem)[i],
                    get<"specific_heating_rate">(chem)[i], get<"RT_heating_rate">(chem)[i],
                    get<"RT_HI_ionization_rate">(chem)[i], get<"RT_HeI_ionization_rate">(chem)[i],
                    get<"RT_HeII_ionization_rate">(chem)[i], get<"RT_H2_dissociation_rate">(chem)[i],
                    get<"H2_self_shielding_length">(chem)[i]);
                d.p[i] = pressure;
            }
        };



        auto calculateChemistry = [&, &chem = simData.chem]()
        {
            const auto oldFields = chem.fields;
            T          max_diff  = 0.;

//#pragma omp parallel for schedule(static) reduction(max : max_diff)
            for (size_t i = first; i < last; ++i)
            {
                T u = d.u[i];
                std::cout << "rho: " << d.rho[i] << "u " << u << std::endl;
                cooling_data.cool_particle(
                    0.01, d.rho[i], u, get<"HI_fraction">(chem)[i], get<"HII_fraction">(chem)[i],
                    get<"HM_fraction">(chem)[i], get<"HeI_fraction">(chem)[i], get<"HeII_fraction">(chem)[i],
                    get<"HeIII_fraction">(chem)[i], get<"H2I_fraction">(chem)[i], get<"H2II_fraction">(chem)[i],
                    get<"DI_fraction">(chem)[i], get<"DII_fraction">(chem)[i], get<"HDI_fraction">(chem)[i],
                    get<"e_fraction">(chem)[i], get<"metal_fraction">(chem)[i], get<"volumetric_heating_rate">(chem)[i],
                    get<"specific_heating_rate">(chem)[i], get<"RT_heating_rate">(chem)[i],
                    get<"RT_HI_ionization_rate">(chem)[i], get<"RT_HeI_ionization_rate">(chem)[i],
                    get<"RT_HeII_ionization_rate">(chem)[i], get<"RT_H2_dissociation_rate">(chem)[i],
                    get<"H2_self_shielding_length">(chem)[i]);
                for (size_t j = 0; j < chem.fields.size(); j++)
                {
                    const T diff = std::abs(chem.fields[j][i] - oldFields[j][i]);
                    max_diff     = std::max(max_diff, diff);
                }
            }
            return max_diff;
        };

        size_t n_it = 0;
        while (true)
        {
            calculatePressure();
            bool good =  (initDependent(simData));
            //const T max_diff = calculateChemistry();
            const T max_diff=0.;
            n_it++;
            std::cout << "equilibrated " << n_it << "\t" << max_diff << std::endl;
            if (max_diff < 1e-6 && good) break;
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

        contractRhoProfileCloud(d.x, d.y, d.z);
        syncCoords<KeyType>(rank, numRanks, numParticlesGlobal, d.x, d.y, d.z, globalBox);

        initPressure(simData, rank, numRanks, globalBox, numParticlesGlobal);

        return globalBox;
    }

    const std::map<std::string, double>& constants() const override { return settings_; }
};

} // namespace sphexa
