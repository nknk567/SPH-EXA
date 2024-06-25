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

template<typename Tpos, typename Tu, typename Ts, typename Tdu, typename Trho>
__global__ void betaCoolingGPUKernel(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Tdu* du,
                                     const Tu* u, Ts star_mass, Ts star_pos_x, Ts star_pos_y, Ts star_pos_z, Ts beta,
                                     Tpos g, const Trho* rho, Ts u_floor, Trho cooling_rho_limit)

{
    cstone::LocalIndex i = first + blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= last) { return; }
    if (rho[i] >= cooling_rho_limit || u[i] <= u_floor) return;

    const double dx    = x[i] - star_pos_x;
    const double dy    = y[i] - star_pos_y;
    const double dz    = z[i] - star_pos_z;
    const double dist2 = dx * dx + dy * dy + dz * dz;
    const double dist  = sqrt(dist2);
    const double omega = sqrt(g * star_mass / (dist2 * dist));
    du[i] += -u[i] * omega / beta;
}

template<typename Tpos, typename Tu, typename Ts, typename Tdu, typename Trho>
void betaCoolingGPU(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, const Tu* u, Tdu* du,
                    Ts star_mass, const Ts* star_pos, Ts beta, Tpos g, const Trho* rho, Ts u_floor,
                    Trho cooling_rho_limit)
{
    cstone::LocalIndex numParticles = last - first;
    unsigned           numThreads   = 256;
    unsigned           numBlocks    = (numParticles + numThreads - 1) / numThreads;

    betaCoolingGPUKernel<<<numBlocks, numThreads>>>(first, last, x, y, z, du, u, star_mass, star_pos[0], star_pos[1],
                                                    star_pos[2], beta, g, rho, u_floor, cooling_rho_limit);

    checkGpuErrors(cudaDeviceSynchronize());
}

template void betaCoolingGPU(size_t, size_t, const double*, const double*, const double* z, const double*, double*, double,
                             const double*, double, double, const float*, double, float);

template<typename Tu, typename Tdu>
struct f
{
    __device__ double operator()(Tu u, Tdu du) { return abs(0.25 * u / du); }
};

template<typename Tu, typename Tdu>
double computeHeatingTimestepGPU(size_t first, size_t last, const Tu* u, const Tdu* du)
{
    return 0.;
    /*  double minDt = 0.25 * thrust::transform_reduce(u + first, u + last, du + first,
                                                     std::numeric_limits<double>::infinity(), f<Tu, Tdu>{},
                                                     thrust::minimum<double>());

      checkGpuErrors(cudaDeviceSynchronize());
      return minDt;*/
}

template double computeHeatingTimestepGPU(size_t, size_t, const double*, const double*);
