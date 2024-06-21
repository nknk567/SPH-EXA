//
// Created by Noah Kubli on 17.04.2024.
//

#pragma once
#include <cmath>
#include "cstone/tree/accel_switch.hpp"
#include "betaCooling_gpu.hpp"

namespace planet
{

template<typename Tpos, typename Tu, typename Ts, typename Tdu, typename Trho>
void betaCoolingImpl(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Tu* u, Tdu* du,
                     Ts star_mass, const Ts* star_pos, Ts beta, Tpos g, Trho* rho, Trho cooling_rho_limit = 1.683e-3)
{
    for (size_t i = first; i < last; i++)
    {
        if (u[i] < 2.5e-5)
        {
            u[i] = 2.5e-5;
            du[i] = std::max(du[i], 0.);
            continue;
        }
        if (rho[i] > cooling_rho_limit) continue;
        const double dx    = x[i] - star_pos[0];
        const double dy    = y[i] - star_pos[1];
        const double dz    = z[i] - star_pos[2];
        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double dist  = std::sqrt(dist2);
        const double omega = std::sqrt(g * star_mass / (dist2 * dist));
        du[i] += -u[i] * omega / beta;
    }
}

template<typename Dataset, typename StarData>
void betaCooling(Dataset& d, size_t startIndex, size_t endIndex, const StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        //betaCoolingGPU(startIndex, endIndex, rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z),
        //               rawPtr(d.devData.u), rawPtr(d.devData.du), star.m, star.position.data(), star.beta, d.g,
          //             rawPtr(d.devData.rho), star.cooling_rho_limit);
        transferToHost(d, startIndex, endIndex, {"rho", "u", "du"});
        betaCoolingImpl(startIndex, endIndex, d.x.data(), d.y.data(), d.z.data(), d.u.data(), d.du.data(), star.m,
                        star.position.data(), star.beta, d.g, d.rho.data(), star.cooling_rho_limit);
        transferToDevice(d, startIndex, endIndex, {"u", "du"});
    }
    else
    {
        betaCoolingImpl(startIndex, endIndex, d.x.data(), d.y.data(), d.z.data(), d.u.data(), d.du.data(), star.m,
                        star.position.data(), star.beta, d.g, d.rho.data(), star.cooling_rho_limit);
    }
}
} // namespace planet
