//
// Created by Noah Kubli on 11.03.2024.
//
#include <cub/cub.cuh>
#include <thrust/execution_policy.h>
#include <thrust/reduce.h>

#include "cstone/cuda/cuda_utils.cuh"
#include "cstone/findneighbors.hpp"
#include "cstone/traversal/find_neighbors.cuh"
#include "sph/util/device_math.cuh"

#include "cstone/sfc/box.hpp"

#include "computeCentralForce_gpu.hpp"
#include "cuda_runtime.h"

template<size_t numThreads, typename Tpos, typename Ta, typename Tm, typename Ts>
__global__ void computeCentralForceGPUKernel(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z,
                                             Ta* ax, Ta* ay, Ta* az, const Tm* m, Ts star_pos_x, Ts star_pos_y,
                                             Ts star_pos_z, Ts sm, Tpos g, Ts* star_force_block_x,
                                             Ts* star_force_block_y, Ts* star_force_block_z, Ts* star_potential_block)
{
    __shared__ Ts star_force_thread_x[numThreads];
    __shared__ Ts star_force_thread_y[numThreads];
    __shared__ Ts star_force_thread_z[numThreads];
    __shared__ Ts star_potential_thread[numThreads];

    cstone::LocalIndex i = first + blockDim.x * blockIdx.x + threadIdx.x;

    if (i >= last)
    {
        star_force_thread_x[threadIdx.x]   = 0.;
        star_force_thread_y[threadIdx.x]   = 0.;
        star_force_thread_z[threadIdx.x]   = 0.;
        star_potential_thread[threadIdx.x] = 0.;
    }
    else
    {
        const double dx    = x[i] - star_pos_x;
        const double dy    = y[i] - star_pos_y;
        const double dz    = z[i] - star_pos_z;
        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double dist  = sqrt(dist2);
        const double dist3 = dist2 * dist;

        const double a_strength = 1. / dist3 * sm * g;
        const double ax_i       = -dx * a_strength;
        const double ay_i       = -dy * a_strength;
        const double az_i       = -dz * a_strength;
        ax[i] += ax_i;
        ay[i] += ay_i;
        az[i] += az_i;

        star_force_thread_x[threadIdx.x]   = -ax_i * m[i];
        star_force_thread_y[threadIdx.x]   = -ay_i * m[i];
        star_force_thread_z[threadIdx.x]   = -az_i * m[i];
        star_potential_thread[threadIdx.x] = -g * m[i] / dist;
    }

    __syncthreads();

    for (unsigned s = blockDim.x / 2; s > 0; s /= 2)
    {
        if (threadIdx.x < s)
        {
            star_force_thread_x[threadIdx.x] += star_force_thread_x[threadIdx.x + s];
            star_force_thread_y[threadIdx.x] += star_force_thread_y[threadIdx.x + s];
            star_force_thread_z[threadIdx.x] += star_force_thread_z[threadIdx.x + s];
            star_potential_thread[threadIdx.x] += star_potential_thread[threadIdx.x + s];
        }
        __syncthreads();
    }

    if (threadIdx.x == 0)
    {
        star_force_block_x[blockIdx.x]   = star_force_thread_x[0];
        star_force_block_y[blockIdx.x]   = star_force_thread_y[0];
        star_force_block_z[blockIdx.x]   = star_force_thread_z[0];
        star_potential_block[blockIdx.x] = star_potential_thread[0];
    }
}

template<typename Tpos, typename Ta, typename Tm, typename Ts>
void computeCentralForceGPU(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Ta* ax, Ta* ay,
                            Ta* az, const Tm* m, const Ts* star_pos, Ts star_mass, Ts* star_force_local,
                            Ts* star_pot_local, Tpos g)
{
    cstone::LocalIndex numParticles = last - first;
    constexpr unsigned numThreads   = 256;
    unsigned           numBlocks    = (numParticles + numThreads - 1) / numThreads;

    star_force_local[0] = 0.;
    star_force_local[1] = 0.;
    star_force_local[2] = 0.;
    *star_pot_local     = 0.;

    Ts* star_force_block_x;
    cudaMalloc(&star_force_block_x, sizeof(Ts) * numBlocks);
    Ts* star_force_block_y;
    cudaMalloc(&star_force_block_y, sizeof(Ts) * numBlocks);
    Ts* star_force_block_z;
    cudaMalloc(&star_force_block_z, sizeof(Ts) * numBlocks);
    Ts* star_pot_block;
    cudaMalloc(&star_pot_block, sizeof(Ts) * numBlocks);

    printf("numParticles: %u\t last: %u\t first: %u\t numBlocks: %u\t numThreads:%u\n", numParticles, last, first,
           numBlocks, numThreads);
    computeCentralForceGPUKernel<numThreads><<<numBlocks, numThreads>>>(
        first, last, x, y, z, ax, ay, az, m, star_pos[0], star_pos[1], star_pos[2], star_mass, g, star_force_block_x,
        star_force_block_y, star_force_block_z, star_pot_block);
    checkGpuErrors(cudaGetLastError());
    checkGpuErrors(cudaDeviceSynchronize());

    star_force_local[0] =
        thrust::reduce(thrust::device, star_force_block_x, star_force_block_x + numBlocks, 0., thrust::plus<Ts>{});
    star_force_local[1] =
        thrust::reduce(thrust::device, star_force_block_y, star_force_block_y + numBlocks, 0., thrust::plus<Ts>{});
    star_force_local[2] =
        thrust::reduce(thrust::device, star_force_block_z, star_force_block_z + numBlocks, 0., thrust::plus<Ts>{});
    *star_pot_local =
        thrust::reduce(thrust::device, star_pot_block, star_pot_block + numBlocks, 0., thrust::plus<Ts>{});

    cudaFree(star_force_block_x);
    cudaFree(star_force_block_y);
    cudaFree(star_force_block_z);
    cudaFree(star_pot_block);
    checkGpuErrors(cudaDeviceSynchronize());
}

template void computeCentralForceGPU(size_t, size_t, const double*, const double*, const double*, float*, float*,
                                     float*, const float*, const double*, double, double*, double*, double);
template void computeCentralForceGPU(size_t, size_t, const double*, const double*, const double*, double*, double*,
                                     double*, const double*, const double*, double, double*, double*, double);
