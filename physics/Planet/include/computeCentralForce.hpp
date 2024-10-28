//
// Created by Noah Kubli on 04.03.2024.
//

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <mpi.h>
#include "cstone/primitives/mpi_wrappers.hpp"
#include "cstone/fields/field_get.hpp"
#include "cstone/tree/accel_switch.hpp"
#include "cstone/cuda/cuda_stubs.h"

#include "computeCentralForce_gpu.hpp"
#include "accretion_gpu.hpp"
#include "accretion_impl.hpp"
#include "sph/particles_data.hpp"

namespace planet
{

template<typename Tpos, typename Ta, typename Tm, typename Ts>
void computeCentralForceImpl(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Ta* ax, Ta* ay,
                             Ta* az, const Tm* m, const Ts* star_pos, Ts star_mass, Ts* star_force_local,
                             Ts* star_pot_local, Tpos g, Ts inner_size)
{
    star_force_local[0]      = 0.;
    star_force_local[1]      = 0.;
    star_force_local[2]      = 0.;
    *star_pot_local          = 0.;
    const double inner_size2 = inner_size * inner_size;

#pragma omp parallel for reduction(+ : star_force_local[ : 3]) reduction(+ : star_pot_local[ : 1])
    for (size_t i = first; i < last; i++)
    {
        const double dx    = x[i] - star_pos[0];
        const double dy    = y[i] - star_pos[1];
        const double dz    = z[i] - star_pos[2];
        const double dist2 = std::max(inner_size2, dx * dx + dy * dy + dz * dz);
        const double dist  = std::sqrt(dist2);
        const double dist3 = dist2 * dist;

        const double a_strength = 1. / dist3 * star_mass * g;
        const double ax_i       = -dx * a_strength;
        const double ay_i       = -dy * a_strength;
        const double az_i       = -dz * a_strength;
        ax[i] += ax_i;
        ay[i] += ay_i;
        az[i] += az_i;

        star_force_local[0] -= ax_i * m[i];
        star_force_local[1] -= ay_i * m[i];
        star_force_local[2] -= az_i * m[i];
        *star_pot_local -= g * m[i] / dist;
    }
}

template<typename Dataset, typename StarData>
void computeCentralForce(Dataset& d, size_t startIndex, size_t endIndex, StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        computeCentralForceGPU(startIndex, endIndex, rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z),
                               rawPtr(d.devData.ax), rawPtr(d.devData.ay), rawPtr(d.devData.az), rawPtr(d.devData.m),
                               star.position.data(), star.m, star.force_local.data(), &star.potential_local, d.g,
                               star.inner_size);
    }
    else
    {
        computeCentralForceImpl(startIndex, endIndex, d.x.data(), d.y.data(), d.z.data(), d.ax.data(), d.ay.data(),
                                d.az.data(), d.m.data(), star.position.data(), star.m, star.force_local.data(),
                                &star.potential_local, d.g, star.inner_size);
    }
}

} // namespace planet
