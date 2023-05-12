//
// Created by Noah Kubli on 04.04.23.
//

#ifndef SPHEXA_STAR_FORMATION_HPP
#define SPHEXA_STAR_FORMATION_HPP

namespace star_formation {

template <typename T, typename Tmass>
struct StarFormer {
    const T temp_max{15000.};
    const T rho_overdensity{20.};
    const T rho_min{0.1};
    const T star_efficiency{0.3333};
    const T t_starform_min{1e6};
    const T c_star{0.05};
    const T m_star{0.0};

    void _createStarParticle(const T& time, const Tmass& mass, const T& metallicity) {
        //Create a new particle and make its type a star
        //...
        //Set new parameters
        //get<"m">(star) = mass;
        //get<"m_form">(star) = mass;
        //get<"t_form">(star) = time;
        //get<"metallicity">(star) = metallicity;
    };

    void form(const T& rho, const T& temp, Tmass &m, const T& dt, const T& time, const T& metallicity)
    {
        constexpr T G{1.0};
        const T t_dynamical = 1. / std::sqrt(4. * M_PI * G * rho);
        if (temp > temp_max) return;
        if (rho < rho_overdensity) return;
        if (rho < rho_min) return;
        const T delta_t = std::max(dt, t_starform_min);

        const float p = 1.0 - std::exp(c_star * delta_t / t_dynamical);

        const Tmass m_star = std::min(m * static_cast<Tmass>(star_efficiency), m);

        const double r_random = std::rand() / (double)RAND_MAX;
        if (m / m_star < r_random / p) return;
        _createStarParticle(time, m_star, metallicity);
        m -= m_star;
        //if (m < m_gas_min) removeParticle();

    }

};
template <typename Dataset>
void form_all(size_t first, size_t last, Dataset &d) {
    using T = typename Dataset::HydroData::RealType;
    using Tmass = typename Dataset::HydroData::Tmass;

    StarFormer<T, Tmass> s{};
    for (size_t i = first; i < last; i++) {
        s.form(get<"rho">(d.hydro)[i],
            get<"temp">(d.hydro)[i],
                get<"m">(d.hydro)[i],
            d.hydro.minDt,
               d.hydro.ttot,
               get<"metal_fraction">(d.chem)[i]);
    }
}

}

#endif // SPHEXA_STAR_FORMATION_HPP
