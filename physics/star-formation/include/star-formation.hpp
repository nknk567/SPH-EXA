//
// Created by Noah Kubli on 04.04.23.
//

#ifndef SPHEXA_STAR_FORMATION_HPP
#define SPHEXA_STAR_FORMATION_HPP

#include <optional>
#include <random>

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
    m -= star_mass;
    if (m < params.m_gas_min) deleteParticle(i);
}

template<typename Dataset>
void form_all(size_t first, size_t last, Dataset& d)
{
    using T     = typename Dataset::HydroData::RealType;
    using Tmass = typename Dataset::HydroData::Tmass;
    const Params<T, Tmass> params;

    for (size_t i = first; i < last; i++)
    {
        form_new(params, i, get<"m">(d.hydro)[i], get<"temp">(d.hydro)[i], get<"divv">(d.hydro)[i],
                 get<"rho">(d.hydro)[i], 1., get<"metal_fraction">(d.chem)[i], get<"oxygen_fraction">(d.chem)[i],
                 get<"iron_fraction">(d.chem)[i], d.hydro.minDt, d.hydro.ttot);
    }
    // redistributeGas();
    /*
     * sum kernel contributions to all neighbours (non-deleted)
     * distribute mass according to fraction of kernel
     * (sph.cpp: 1993)
     */
}

} // namespace star_formation

#endif // SPHEXA_STAR_FORMATION_HPP
