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

template<typename Dataset, typename StarData>
void betaCoolingImpl(size_t first, size_t last, Dataset& d, const StarData& star)
{
#pragma omp parallel for
    for (size_t i = first; i < last; i++)
    {
        if (d.rho[i] < star.cooling_rho_limit && d.u[i] > star.u_floor)
        {
            const double dx    = d.x[i] - star.position[0];
            const double dy    = d.y[i] - star.position[1];
            const double dz    = d.z[i] - star.position[2];
            const double dist2 = dx * dx + dy * dy + dz * dz;
            const double dist  = std::sqrt(dist2);
            const double omega = std::sqrt(d.g * star.m / (dist2 * dist));
            d.du[i] += -d.u[i] * omega / star.beta;
        }
    }
}

// template<std::floating_point Tt>
// Tt duTimestepAndTempFloorImpl(size_t first, size_t last, std::floating_point auto* du, std::floating_point auto* u,
//                               const std::floating_point auto* du_m1, std::floating_point auto u_floor,
//                               std::floating_point auto u_max, Tt k_u)
//{
//     size_t n_below_floor{};
//     size_t n_above_max{};
//     Tt     duTimestepMin = std::numeric_limits<Tt>::infinity();
//
// #pragma omp parallel for reduction(min : duTimestepMin) reduction(+ : n_below_floor) reduction(+ : n_above_max)
//     for (size_t i = first; i < last; i++)
//     {
//         if (u[i] < u_floor)
//         {
//             u[i]  = u_floor;
//             du[i] = std::max(0., du[i]);
//             n_below_floor++;
//         }
//         else if (u[i] > u_max)
//         {
//             u[i]  = u_max;
//             du[i] = std::min(0., du[i]);
//             n_above_max++;
//         }
//
//         Tt duTimestep = k_u * std::abs(u[i] / du[i]);
//         duTimestepMin = std::min(duTimestepMin, duTimestep);
//     }
//     printf("n_below_floor: %zu\t n_above_max: %zu\n", n_below_floor, n_above_max);
//     return duTimestepMin;
// }

template<typename Dataset, typename StarData>
auto duTimestepImpl(size_t first, size_t last, const Dataset& d, const StarData& star)
{

    using Tu         = std::decay_t<decltype(d.u[0])>;
    using Tdu        = std::decay_t<decltype(d.du[0])>;
    using Tt         = std::common_type_t<Tu, Tdu>;
    Tt duTimestepMin = std::numeric_limits<Tt>::infinity();

#pragma omp parallel for reduction(min : duTimestepMin)
    for (size_t i = first; i < last; i++)
    {
        Tt duTimestep = star.K_u * std::abs(d.u[i] / d.du[i]);
        duTimestepMin = std::min(duTimestepMin, duTimestep);
    }
    return duTimestepMin;
}

template<typename Dataset, typename StarData>
void betaCooling(size_t startIndex, size_t endIndex, Dataset& d, const StarData& star)
{
    using T_beta = std::decay_t<decltype(star.beta)>;
    if (star.beta != std::numeric_limits<T_beta>::infinity())
    {
        if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
        {
            betaCoolingGPU(startIndex, endIndex, rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z),
                           rawPtr(d.devData.u), rawPtr(d.devData.du), star.m, star.position.data(), star.beta, d.g,
                           rawPtr(d.devData.rho), star.u_floor, star.cooling_rho_limit);
        }
        else { betaCoolingImpl(startIndex, endIndex, d, star); }
    }
}

template<typename Dataset, typename StarData>
void duTimestep(size_t startIndex, size_t endIndex, Dataset& d, StarData& star)
{
    if (star.u_floor == 0. && star.u_max == std::numeric_limits<decltype(star.u_max)>::infinity() &&
        star.du_adjust == std::numeric_limits<decltype(star.du_adjust)>::infinity())
    {
        return;
    }
    else
    {
        if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
        {
            star.t_du = duTimestepGPU(startIndex, endIndex, d, star);
        }
        else { star.t_du = duTimestepImpl(startIndex, endIndex, d, star); }
    }
}
// template<typename Dataset, typename StarData>
// void duTimestepAndTempFloor(Dataset& d, size_t startIndex, size_t endIndex, StarData& star)
//{
//     if (star.u_floor == 0. && star.u_max == std::numeric_limits<decltype(star.u_max)>::infinity() &&
//         star.du_adjust == std::numeric_limits<decltype(star.du_adjust)>::infinity())
//     {
//         return;
//     }
//     else
//     {
//         if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
//         {
//             transferToHost(d, startIndex, endIndex, {"du", "u"});
//
//             auto dt_u = duTimestepAndTempFloorImpl(startIndex, endIndex, d.du.data(), d.u.data(), d.du_m1.data(),
//                                                    star.u_floor, star.u_max, star.K_u);
//
//             transferToDevice(d, startIndex, endIndex, {"du", "u"});
//             star.t_du = dt_u;
//         }
//         else
//         {
//             auto dt_u = duTimestepAndTempFloorImpl(startIndex, endIndex, d.du.data(), d.u.data(), d.du_m1.data(),
//                                                    star.u_floor, star.u_max, star.K_u);
//             star.t_du = dt_u;
//         }
//     }
// }
} // namespace planet
