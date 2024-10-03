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

#include "cstone/domain/domain.hpp"
#include "cstone/sfc/box.hpp"
#include "cstone/tree/continuum.hpp"
#include "cstone/util/constexpr_string.hpp"
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
    return {{"gravConstant", 1.},
            {"r", 1.5},
            {"mTotal", 1.},
            {"gamma", 5. / 3.},
            {"u0", 0.05},
            {"minDt", 1e-4},
            {"minDt_m1", 1e-4},
            {"mui", 10},
            {"ng0", 100},
            {"ngmax", 110},
            {"metal_fraction", 1e-6},
            {"hydrogen_fraction", 0.76},
            {"d_to_h_ratio", 1e-5},
            {"cooling::ct_crit", 0.01},
            {"cooling::use_grackle", 1},
            {"cooling::with_radiative_cooling", 1},
            {"cooling::primordial_chemistry", 1},
            {"cooling::dust_chemistry", 0},
            {"cooling::metal_cooling", 0},
            {"cooling::UVbackground", 0},
            {"cooling::m_code_in_ms", 1e8},
            {"cooling::l_code_in_kpc", 1.0},
            {"cooling::collisional_excitation_rates", 1},
            {"cooling::collisional_ionisation_rates", 1},
            {"cooling::recombination_cooling_rates", 1},
            {"cooling::bremsstrahlung_cooling_rates", 1},
            {"cooling::use_specific_heating_rate", 1}};
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

    // std::fill(d.soft.begin(), d.soft.end(), 0.05);

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
        T    radius        = std::sqrt((d.x[i] * d.x[i]) + (d.y[i] * d.y[i]) + (d.z[i] * d.z[i]));
        auto concentration = [c0](double x, double y, double z)
        {
            const double r = std::sqrt(x * x + y * y + z * z);
            if (r < 1.)
                return c0;
            else
                return c0 / std::pow(r, 4);
        };
        d.h[i] = std::cbrt(3 / (4 * M_PI) * d.ng0 / concentration(d.x[i], d.y[i], d.z[i])) * 0.5;
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

//#pragma omp parallel for schedule(static)
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
    mutable InitSettings settings_;

public:
    explicit CloudGlassSphere(std::string initBlock, std::string settingsFile, IFileReader* reader)
        : glassBlock(std::move(initBlock))
    {
        Dataset d;
        settings_ = buildSettings(d, cloudConstants(), settingsFile, reader);
        resetConstants(settings_);
        // Base::updateSettings(cloudConstants());
    }
    void resetConstants(InitSettings newSettings) { settings_ = std::move(newSettings); }

    bool adjustInternalEnergy(Dataset& simData) const
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

    template<typename Domain, typename Cooler>
    void calculatePressure(Dataset& simData, Domain& domain, Cooler& cooling_data) const
    {
        using T            = typename Dataset::RealType;
        auto&        d     = simData.hydro;
        auto&        chem  = simData.chem;
        const size_t first = domain.startIndex();
        const size_t last  = domain.endIndex();
        using ChemData     = typename Dataset::ChemData;

        using CoolingFields = typename cooling::Cooler<T>::CoolingFields;
        cooling_data.computePressures(d.rho.data(), d.u.data(), cstone::getPointers(get<CoolingFields>(chem), 0),
                                      d.p.data(), first, last);
    };

    template<typename Domain, typename Cooler>
    void calculatePressureIdealMonatomic(Dataset& simData, Domain& domain, Cooler& cooling_data) const
    {
        using T            = typename Dataset::RealType;
        auto&        d     = simData.hydro;
        const size_t first = domain.startIndex();
        const size_t last  = domain.endIndex();
        constexpr T  gamma = 5. / 3.;
#pragma omp parallel for schedule(static)
        for (size_t i = first; i < last; ++i)
        {
            d.p[i] = (gamma - 1.) * d.u[i] * d.rho[i];
        }
    };

    template<typename Domain, typename Cooler>
    typename Dataset::RealType calculateChemistry(Dataset& simData, Domain& domain, Cooler& cooling_data,
                                                  typename Dataset::RealType timestep) const
    {
        using T            = typename Dataset::RealType;
        auto&        d     = simData.hydro;
        auto&        chem  = simData.chem;
        const size_t first = domain.startIndex();
        const size_t last  = domain.endIndex();
        using ChemData     = typename Dataset::ChemData;
        //        using CoolingFields = typename util::MakeFieldList<ChemData>::Fields;
        using CoolingFields = typename cooling::Cooler<T>::CoolingFields;

        T max_diff = 0.;

#pragma omp parallel for schedule(static) reduction(max : max_diff)
        for (size_t i = first; i < last; ++i)
        {
            std::array<T, chem.numFields> old;
            for (size_t j = 0; j < chem.fields.size(); j++)
            {
                old[j] = chem.fields[j][i];
            }
            T u   = d.u[i];
            T rho = d.rho[i];
            cooling_data.cool_particles(timestep, d.rho.data(), d.u.data(),
                                        cstone::getPointers(get<CoolingFields>(chem), 0), d.du.data(), i, i + 1);
            for (size_t j = 0; j < chem.fields.size(); j++)
            {
                const T diff = std::abs(chem.fields[j][i] - old[j]) / timestep;
                max_diff     = std::max(max_diff, diff);
            }
        }
        return max_diff;
    };

    void initPressure(Dataset& simData, const int rank, const int numRanks,
                      const cstone::Box<typename Dataset::RealType>& globalBox, const size_t numParticlesGlobal) const
    {
        auto& d       = simData.hydro;
        using KeyType = typename Dataset::KeyType;
        using T       = typename Dataset::RealType;
        // using CoolingFields = typename util::MakeFieldList<ChemData>::Fields;
        using util::FieldList;
        using ConservedFields = FieldList<"temp", "vx", "vy", "vz", "x_m1", "y_m1", "z_m1", "du_m1", "u">;

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
        auto settings_init                               = settings_;
        settings_init["cooling::with_radiative_cooling"] = 0;
        settings_init["cooling::primordial_chemistry"]   = 1;
        settings_init["cooling::max_iterations"]         = 10000 * 10000;
        BuiltinWriter attributeSetter(settings_init);
        d.loadOrStoreAttributes(&attributeSetter);
        cooling::Cooler<T> cooling_data;
        cooling_data.loadOrStoreAttributes(&attributeSetter);

        initCloudFields(d, simData.chem, settings_, globalBox);

        // Prepare hydrost. eq
        uint64_t bucketSizeFocus = 64;

        uint64_t bucketSize = std::max(bucketSizeFocus, d.numParticlesGlobal / (100 * numRanks));
        cstone::Domain<KeyType, T, cstone::CpuTag> domain(rank, numRanks, bucketSize, bucketSizeFocus, 0.5, globalBox);

        if (rank == 0) std::cout << "nLocalParticles " << get<"HI_fraction">(simData.chem).size() << std::endl;

        /*auto bl = []()
        {
            std::vector<std::string> ret{"x", "y", "z", "h", "m"};
            for_each_tuple([&ret](auto f) { ret.push_back(f.value); }, make_tuple(ConservedFields{}));
            return ret;
        };*/

        // transferToDevice(d, 0, d.x.size(), bl());

        domain.syncGrav(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                        std::tuple_cat(get<ConservedFields>(d), get<CoolingFields>(simData.chem)),
                        get<DependentFields>(d));
        d.treeView = domain.octreeProperties();

        d.resize(domain.nParticlesWithHalos());

        resizeNeighbors(d, domain.nParticles() * d.ngmax);
        const size_t first = domain.startIndex();
        const size_t last  = domain.endIndex();

        sph::GroupData<cstone::CpuTag> groups;

        sph::findNeighborsSfc(first, last, d, domain.box());
        sph::computeGroups(first, last, d, domain.box(), groups);

        sph::computeDensity(groups.view(), d, domain.box()); // halo exchange rho!!

        cooling_data.init(0);

        const bool use_grackle = false;
        if (use_grackle)
        {
            size_t n_it     = 0;
            T      timestep = 1.;
            T      time     = 0.;
            while (true)
            {
                calculatePressure(simData, domain, cooling_data);
                bool good = (adjustInternalEnergy(simData));
                auto diff = calculateChemistry(simData, domain, cooling_data, timestep);
                n_it++;
                time += timestep;
                std::cout << "equilibrated " << n_it << std::endl;
                std::cout << "timestep " << timestep << std::endl;
                std::cout << "diff " << diff << std::endl;
                std::cout << "time " << time << std::endl;
                if (good && diff < 0.039 /*1e-4*/) break;
                // if (diff * timestep < 1e-3) timestep *= 1.5;
            }
        }
        else
        {
            // Else
            calculatePressureIdealMonatomic(simData, domain, cooling_data); // oder einfach standard eos verwenden
            bool good = adjustInternalEnergy(simData);
            calculatePressureIdealMonatomic(simData, domain, cooling_data);
            good = (adjustInternalEnergy(simData));
            if (!good) throw std::runtime_error("not good");
        }
    }

    cstone::Box<typename Dataset::RealType> init(int rank, int numRanks, size_t cbrtNumPart, Dataset& simData,
                                                 IFileReader* reader) const override
    {
        auto& d       = simData.hydro;
        using KeyType = typename Dataset::KeyType;
        using T       = typename Dataset::RealType;

        std::vector<T> xBlock, yBlock, zBlock;
        sphexa::readTemplateBlock(glassBlock, reader, xBlock, yBlock, zBlock);
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
printf("r: %lf\n", r);
        initPressure(simData, rank, numRanks, globalBox, numParticlesGlobal);

        return globalBox;
    }

    const std::map<std::string, double>& constants() const override { return settings_; }
};

} // namespace sphexa
