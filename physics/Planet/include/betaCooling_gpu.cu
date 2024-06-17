//
// Created by Noah Kubli on 17.04.2024.
//

#include "betaCooling_gpu.hpp"
#include "sph/util/device_math.cuh"
#include "cstone/cuda/cuda_utils.cuh"
#include "cstone/findneighbors.hpp"
#include "cstone/traversal/find_neighbors.cuh"
#include "sph/util/device_math.cuh"

#include "cstone/sfc/box.hpp"

template<typename Tpos, typename Ts, typename Tdu, typename Trho, typename Tu>
__global__ void betaCoolingGPUKernel(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Tdu* du, const Tu *u,
                                     Ts star_mass, Ts star_pos_x, Ts star_pos_y, Ts star_pos_z, Ts beta, Tpos g,
                                     Trho *rho, Trho cooling_rho_limit)

{
    cstone::LocalIndex i = first + blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= last) { return; }
    if (rho[i] > cooling_rho_limit) return;

    const double dx    = x[i] - star_pos_x;
    const double dy    = y[i] - star_pos_y;
    const double dz    = z[i] - star_pos_z;
    const double dist2 = dx * dx + dy * dy + dz * dz;
    const double dist  = sqrt(dist2);
    const double omega = sqrt(g * star_mass / (dist2 * dist));
    du[i] += -u[i] * omega / beta;
}

template<typename Tpos, typename Ts, typename Tdu, typename Trho, typename Tu>
void betaCoolingGPU(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Tdu* du, const Tu* u, Ts star_mass,
                    const Ts* star_pos, Ts beta, Tpos g, Trho *rho, Trho cooling_rho_limit)
{
    cstone::LocalIndex numParticles = last - first;
    unsigned           numThreads   = 256;
    unsigned           numBlocks    = (numParticles + numThreads - 1) / numThreads;

    betaCoolingGPUKernel<<<numBlocks, numThreads>>>(first, last, x, y, z, du, u, star_mass, star_pos[0], star_pos[1],
                                                    star_pos[2], beta, g, rho, cooling_rho_limit);

    checkGpuErrors(cudaDeviceSynchronize());
}

template void betaCoolingGPU(size_t, size_t, const double*, const double*, const double* z, double*, const double*, double,
                             const double*, double, double, float *, float);
