//
// Created by Noah Kubli on 12.03.2024.
//
#include <cub/cub.cuh>

#include "cstone/cuda/cuda_utils.cuh"
#include "cstone/findneighbors.hpp"
#include "cstone/traversal/find_neighbors.cuh"
//#include "sph/util/device_math.cuh"

#include "cstone/sfc/box.hpp"
#include "cstone/tree/definitions.h"

#include "sph/particles_data.hpp"

#include "accretion_gpu.hpp"
#include "star_data.hpp"
#include "removalStatistics.hpp"
#include "cuda_runtime.h"

namespace disk
{

__device__ void atomicAddRS(RemovalStatistics* x, const RemovalStatistics& y)
{
    atomicAdd(&(x->mass), y.mass);
    atomicAdd(&(x->momentum[0]), y.momentum[0]);
    atomicAdd(&(x->momentum[1]), y.momentum[1]);
    atomicAdd(&(x->momentum[2]), y.momentum[2]);
    atomicAdd(&(x->count), y.count);
}

template<typename Tkeys>
__device__ void markForRemovalAndAdd(RemovalStatistics& statistics, size_t i, Tkeys* keys, const auto* m,
                                     const auto* vx, const auto* vy, const auto* vz)
{
    keys[i]                = cstone::removeKey<Tkeys>::value;
    statistics.mass        = m[i];
    statistics.momentum[0] = m[i] * vx[i];
    statistics.momentum[1] = m[i] * vy[i];
    statistics.momentum[2] = m[i] * vz[i];
    statistics.count       = 1;
}

template<unsigned numThreads, typename T1, typename Th, typename Tkeys, typename T2, typename Tm, typename Tv>
__global__ void computeAccretionConditionKernel(size_t first, size_t last, const T1* x, const T1* y, const T1* z,
                                                const Th* h, Tkeys* keys, const Tm* m, const Tv* vx, const Tv* vy,
                                                const Tv* vz, T2 star_x, T2 star_y, T2 star_z, T2 star_size2,
                                                T2 removal_limit_h, RemovalStatistics* device_accreted,
                                                RemovalStatistics* device_removed)
{
    cstone::LocalIndex i = first + blockDim.x * blockIdx.x + threadIdx.x;

    // Accreted particles statistics
    RemovalStatistics accreted{};
    // Removed particles statistics
    RemovalStatistics removed{};

    if (i >= last) {}
    else
    {
        const double dx    = x[i] - star_x;
        const double dy    = y[i] - star_y;
        const double dz    = z[i] - star_z;
        const double dist2 = dx * dx + dy * dy + dz * dz;

        if (dist2 < star_size2) { markForRemovalAndAdd(accreted, i, keys, m, vx, vy, vz); }
        else if (h[i] > removal_limit_h) { markForRemovalAndAdd(removed, i, keys, m, vx, vy, vz); }
    }

    typedef cub::BlockReduce<RemovalStatistics, numThreads> BlockReduce;
    __shared__ typename BlockReduce::TempStorage            temp_storage;

    RemovalStatistics block_accreted = BlockReduce(temp_storage).Sum(accreted);
    __syncthreads();
    if (threadIdx.x == 0) { atomicAddRS(device_accreted, block_accreted); }

    RemovalStatistics block_removed = BlockReduce(temp_storage).Sum(removed);
    __syncthreads();
    if (threadIdx.x == 0) { atomicAddRS(device_removed, block_removed); }
}

template<typename Dataset>
void computeAccretionConditionGPU(size_t first, size_t last, Dataset& d, StarData& star)
{
    cstone::LocalIndex numParticles = last - first;
    constexpr unsigned numThreads   = 256;
    unsigned           numBlocks    = (numParticles + numThreads - 1) / numThreads;

    star.accreted_local = {};
    star.removed_local  = {};

    RemovalStatistics *accreted_device, *removed_device;
    checkGpuErrors(cudaMalloc(reinterpret_cast<void**>(&accreted_device), sizeof *accreted_device));
    checkGpuErrors(cudaMalloc(reinterpret_cast<void**>(&removed_device), sizeof *removed_device));
    checkGpuErrors(
        cudaMemcpy(accreted_device, &star.accreted_local, sizeof star.accreted_local, cudaMemcpyHostToDevice));
    checkGpuErrors(cudaMemcpy(removed_device, &star.removed_local, sizeof star.removed_local, cudaMemcpyHostToDevice));

    computeAccretionConditionKernel<numThreads><<<numBlocks, numThreads>>>(
        first, last, rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z), rawPtr(d.devData.h),
        rawPtr(d.devData.keys), rawPtr(d.devData.m), rawPtr(d.devData.vx), rawPtr(d.devData.vy), rawPtr(d.devData.vz),
        star.position[0], star.position[1], star.position[2], star.inner_size * star.inner_size, star.removal_limit_h,
        accreted_device, removed_device);

    checkGpuErrors(cudaDeviceSynchronize());
    checkGpuErrors(cudaGetLastError());

    checkGpuErrors(
        cudaMemcpy(&star.accreted_local, accreted_device, sizeof star.accreted_local, cudaMemcpyDeviceToHost));
    checkGpuErrors(cudaMemcpy(&star.removed_local, removed_device, sizeof star.removed_local, cudaMemcpyDeviceToHost));
    checkGpuErrors(cudaFree(accreted_device));
    checkGpuErrors(cudaFree(removed_device));
}

template void computeAccretionConditionGPU(size_t, size_t, sphexa::ParticlesData<cstone::GpuTag>&, StarData&);
} // namespace disk