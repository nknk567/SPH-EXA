//
// Created by Noah Kubli on 17.04.2024.
//

#pragma once
#include <cmath>
#include "cstone/tree/accel_switch.hpp"
#include "betaCooling_gpu.hpp"
#include "sph/particles_data.hpp"
namespace planet
{

template<typename Tpos, typename Ts, typename Tdu, typename Trho, typename Tu>
void betaCoolingImpl(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Tdu* du, Tu* u,
                     Ts star_mass, const Ts* star_pos, Ts beta, Tpos g, Trho* rho, double* timestep,
                     Trho cooling_rho_limit = 1.683e-3)
{
    *timestep = INFINITY;
    double cooling_floor = 6e-6;
    for (size_t i = first; i < last; i++)
    {
        if (rho[i] > cooling_rho_limit) continue;
        const double dx    = x[i] - star_pos[0];
        const double dy    = y[i] - star_pos[1];
        const double dz    = z[i] - star_pos[2];
        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double dist  = std::sqrt(dist2);
        const double omega = std::sqrt(g * star_mass / (dist2 * dist));
        du[i] += -u[i] * omega / beta;

        if (u[i] <= cooling_floor)
        {
            u[i] = cooling_floor;
            du[i] = std::max(0., du[i]);
        }

        *timestep = std::min(*timestep, std::abs(u[i] / du[i]));
    }
    *timestep *= 0.1;
}

template<typename Dataset, typename StarData>
double betaCooling(Dataset& d, size_t startIndex, size_t endIndex, const StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        transferToHost(d, startIndex, endIndex, {"x", "y", "z", "du", "u", "rho"});
        double timestep = 0.;
        betaCoolingImpl(startIndex, endIndex, d.x.data(), d.y.data(), d.z.data(), d.du.data(), d.u.data(), star.m,
                        star.position.data(), star.beta, d.g, d.rho.data(), &timestep, star.cooling_rho_limit);
        transferToDevice(d, startIndex, endIndex, {"du", "u"});
        return timestep;

        /*betaCoolingGPU(startIndex, endIndex, rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z),
                       rawPtr(d.devData.du), rawPtr(d.devData.u), star.m, star.position.data(), star.beta, d.g,
                       rawPtr(d.devData.rho), star.cooling_rho_limit);*/
    }
    else
    {
        double timestep = 0.;
        betaCoolingImpl(startIndex, endIndex, d.x.data(), d.y.data(), d.z.data(), d.du.data(), d.u.data(), star.m,
                        star.position.data(), star.beta, d.g, d.rho.data(), &timestep, star.cooling_rho_limit);
        return timestep;
    }
}
} // namespace planet
