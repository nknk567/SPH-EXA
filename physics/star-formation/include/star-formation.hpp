//
// Created by Noah Kubli on 04.04.23.
//

#ifndef SPHEXA_STAR_FORMATION_HPP
#define SPHEXA_STAR_FORMATION_HPP

#include <optional>
#include <random>

#include "sph/kernels.hpp"
#include "sph/math.hpp"
#include "sph/tables.hpp"

namespace star_formation
{

template<typename T, typename Tmass>
struct Params
{
    T                    temp_max{15000.};
    T                    rho_overdensity{20.};
    T                    rho_min{0.1};
    T                    star_efficiency{0.3333};
    T                    t_starform_min{1e6};
    T                    c_star{0.05};
    std::optional<Tmass> m_star_init{std::nullopt};
    Tmass                m_gas_min{0.0};
    Tmass                m_star_max_hybrid{0.0};
    T                    soft_min{1.0};
    uint64_t             rung{0};
};
template<typename T, typename Tmass>
size_t addParticle(const Tmass mass, const Tmass formation_mass, const T metallicity, const T oxygen_frac,
                   const T iron_frac, const T formation_time)
{
    return 0;
};
void deleteParticle(size_t i) {}

template<typename T, typename Tmass>
bool conditions(const Params<T, Tmass>& params, const T temp, const T divv, const T rho, const T a)
{
    if (temp > params.temp_max) return false;
    if (divv > 0.) return false;
    if (rho < params.rho_overdensity) return false;
    if (rho / a < params.rho_min) return false;
    return true;
}

template<typename T, typename Tmass>
T formationProbability(const Params<T, Tmass>& params, const T rho, const T a, const T dt)
{
    const T dynamical_time{1. / std::sqrt(4. * M_PI * rho / a)};
    const T star_formation_time{std::max(dt, params.t_starform_min)};
    const T p{1. - std::exp(-params.c_star * star_formation_time / dynamical_time)};
    return p;
}

template<typename T, typename Tmass>
T formationMass(const Params<T, Tmass>& params, Tmass& m, const T temp, const T divv, const T rho, const T a,
                const T dt)
{
    if (!conditions(params, temp, divv, rho, a)) return T{0.0};
    const T     p         = formationProbability(params, rho, a, dt);
    const Tmass star_mass = [&]()
    {
        if (params.m_star_init.has_value())
            return params.m_star_init.value();
        else
            return m * Tmass(params.star_efficiency);
    }();
    const Tmass star_mass_constr = std::min(star_mass, m);
    const T     rand_val         = rand() / T(RAND_MAX);
    if (p * m < rand_val * star_mass_constr) return T{0.0};
    return star_mass_constr;
}

template<typename T>
struct particle_type
{
    inline constexpr static T gas{0.0};
    inline constexpr static T star{1.0};
};

template<typename Dataset, typename T>
void initializeStellarParticle(Dataset& d, const size_t ig, const size_t is, const T stellarMass)
{
    // Copy information to star particle
    d.hydro.m[is]                     = stellarMass;
    get<"metal_fraction">(d.chem)[is] = get<"metal_fraction">(d.chem)[ig];
    d.hydro.pType[is] = particle_type<T>::star;
    /*get<"formation_time">(d.star)[is] = d.hydro.ttot;
    get<"formation_mass">(d.star)[is] = stellarMass;*/

    /*  (const Tmass mass, const Tmass formation_mass, const T metallicity, const T oxygen_frac,
       const T iron_frac, const T formation_time)*/

    /* et<"m">(d.hydro)[i], get<"temp">(d.hydro)[i], get<"divv">(d.hydro)[i],
         get<"rho">(d.hydro)[i], 1., get<"metal_fraction">(d.chem)[i], get<"oxygen_fraction">(d.chem)[i],
         get<"iron_fraction">(d.chem)[i], d.hydro.minDt, d.hydro.ttot*/
}

template<typename Domain, typename Dataset>
void resizeStar(Domain& domain, Dataset& d, size_t n_new)
{
    if (n_new > domain.nParticlesWithHalos() - domain.nParticles())
    {
        size_t new_size = domain.nParticles() + n_new;
        d.hydro.resize(new_size);
        d.chem.resize(new_size);
        std::cout << "new_size: " << new_size << std::endl;
    }
    const size_t end_before = domain.endIndex();
    domain.setEndIndex(end_before + n_new);
}

template<typename Domain, typename Dataset>
void form_all(Domain& domain, size_t first, size_t last, Dataset& d)
{
    using T     = typename Dataset::HydroData::RealType;
    using Tmass = typename Dataset::HydroData::Tmass;
    const Params<T, Tmass> params;

    size_t              n_formation        = 0;
    Tmass               tot_formation_mass = 0.;
    std::vector<size_t> formation_indices{};
    std::vector<Tmass>  formation_masses{};

#pragma omp parallel reduction(+ : n_formation)
    {
        std::vector<size_t> local_formation_indices{};
        std::vector<Tmass>  local_formation_masses{};
#pragma omp for nowait // reduction(+: n_formation) reduction(+: tot_formation_mass)
        for (size_t i = first; i < last; i++)
        {
            Tmass formation_mass = formationMass(params, d.hydro.m[i], d.hydro.temp[i], d.hydro.divv[i], d.hydro.rho[i],
                                                 1.0, d.hydro.minDt);
            if (formation_mass > Tmass{0.}) n_formation++;
            formation_masses.push_back(formation_mass);
            formation_indices.push_back(i);
        }
#pragma omp critical
        {
            formation_indices.insert(formation_indices.end(), local_formation_indices.begin(),
                                     local_formation_indices.end());
            formation_masses.insert(local_formation_masses.end(), local_formation_masses.begin(),
                                    local_formation_masses.end());
        }
    }

    resizeStar(domain, d, n_formation);

#pragma omp parallel for
    for (size_t i = 0; i < formation_indices.size(); i++)
    {
        const size_t index{formation_indices[i]};
        const Tmass  m_star{formation_masses[i]};

        initializeStellarParticle(d, index, domain.endIndex() - i - 1, m_star);

        auto& m_gas = d.hydro.m[index];
        m_gas -= m_star;
        if (m_gas < params.m_gas_min) d.hydro.redistributeParticle[index] = 1.0;
    }
}
// redistributeGas();
/*
 * sum kernel contributions to all neighbours (non-deleted)
 * distribute mass according to fraction of kernel
 * (sph.cpp: 1993)
 */
//}

/*template<typename T, typename Tmass>
void form_new(const Params<T, Tmass>& params, const size_t i, Tmass& m, const T temp, const T divv, const T rho,
              const T a, const T metallicity, const T oxygen_frac, const T iron_frac, const T dt, const T sim_time)
{
    if (!conditions(params, temp, divv, rho, a)) return;
    const T     p         = formationProbability(params, rho, a, dt);
    const Tmass star_mass = [&]()
    {
        if (params.m_star_init.has_value())
            return params.m_star_init.value();
        else
            return m * Tmass(params.star_efficiency);
    }();
    const Tmass star_mass_constr = std::min(star_mass, m);
    const T     rand_val         = rand() / T(RAND_MAX);
    if (p * m < rand_val * star_mass_constr) return;

    addParticle(star_mass_constr, star_mass_constr, metallicity, oxygen_frac, iron_frac, sim_time);
    m -= star_mass; //not star_mass_constr?
    if (m < params.m_gas_min)
        deleteParticle(i); // Mark particle to be deleted and calculate w_tot (sum of kernel weights)
}*/

// template<typename T, typename Tc>
// void redistributeDeletedParticle_reversed(cstone::LocalIndex i, T sincIndex, const cstone::LocalIndex* neighbors,
//                                           unsigned neighborsCount, const Tc* x, const Tc* y, const Tc* z, const T* h,
//                                           const T* m, const T* w_tot, const T* is_deleted, const T* wh, const T* whd)
//{
//     auto  xi    = x[i];
//     auto  yi    = y[i];
//     auto  zi    = z[i];
//     auto& mi    = m[i];
//     auto  hi    = h[i];
//     auto  hiInv = T(1) / hi;
//
//     for (unsigned pj = 0; pj < neighborsCount; ++pj)
//     {
//         cstone::LocalIndex j = neighbors[/*stride * */ pj];
//         if (is_deleted[j])
//         {
//             const T rx      = (xi - x[j]);
//             const T ry      = (yi - y[j]);
//             const T rz      = (zi - z[j]);
//             const T dist    = std::sqrt(rx * rx + ry * ry + rz * rz);
//             const T vloc    = dist * hiInv;
//             T       w       = sph::math::pow(sph::lt::wharmonic_lt_with_derivative(wh, whd, vloc), (int)sincIndex);
//             const T m_delta = m[j] * w / w_tot[j];
//
//             const T m_tot   = mi + m_delta;
//             const T c_delta = m_delta / m_tot;
//             const T ci      = mi / m_tot;
//
//             mi = m_tot;
//         }
//     }
// }
// template<typename T, typename Tc>
// void redistributeDeletedParticle(cstone::LocalIndex i, T sincIndex, const cstone::LocalIndex* neighbors,
//                                  unsigned neighborsCount, const Tc* x, const Tc* y, const Tc* z, const T* h, const T*
//                                  m, const T* wh, const T* whd)
//{
//     auto xi    = x[i];
//     auto yi    = y[i];
//     auto zi    = z[i];
//     auto mi    = m[i];
//     auto hi    = h[i];
//     auto hiInv = T(1) / hi;
//
//     std::vector<Tc> w_neighbours(neighborsCount);
//
//     for (unsigned pj = 0; pj < neighborsCount; ++pj)
//     {
//         cstone::LocalIndex j = neighbors[/*stride * */ pj];
//
//         const T rx = (xi - x[j]);
//         const T ry = (yi - y[j]);
//         const T rz = (zi - z[j]);
//         // applyPBC(box, T(2) * hi, rx, ry, rz);
//         const T dist     = std::sqrt(rx * rx + ry * ry + rz * rz);
//         const T vloc     = dist * hiInv;
//         T       w        = sph::math::pow(sph::lt::wharmonic_lt_with_derivative(wh, whd, vloc), (int)sincIndex);
//         w_neighbours[pj] = w;
//     }
//     const T w_tot{std::accumulate(w_neighbours.begin(), w_neighbours.end(), Tc{0.})};
//     for (unsigned pj = 0; pj < neighborsCount; ++pj)
//     {
//         cstone::LocalIndex j = neighbors[pj];
//
//         const T m_delta = mi * w_neighbours[pj] / w_tot;
//
//         const T m_tot   = m[j] + m_delta;
//         const T c_delta = m_delta / m_tot;
//         const T cj      = m[j] / m_tot;
//
//         m[pj] = m_tot;
//         // Velocity, internal energy / Temperature, mass fractions
//     }
// }
//
// template<typename Dataset>
// void redistributeGas(Dataset& d, std::vector<size_t> particles)
//{
// }

} // namespace star_formation

#endif // SPHEXA_STAR_FORMATION_HPP
