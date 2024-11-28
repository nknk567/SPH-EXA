//
// Created by Noah Kubli on 12.03.2024.
//
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/partition.h>
#include <thrust/sequence.h>

#include "cstone/cuda/cuda_utils.cuh"
#include "cstone/findneighbors.hpp"
#include "cstone/traversal/find_neighbors.cuh"
#include "sph/util/device_math.cuh"

#include "cstone/sfc/box.hpp"
#include "cstone/tree/definitions.h"

#include "sph/particles_data.hpp"

#include "accretion_gpu.hpp"
#include "star_data.hpp"
#include "cuda_runtime.h"

struct RemovalProperties
{
    double   mass        = 0.;
    double   momentum[3] = {0.};
    unsigned count       = 0;
    //    RemovalProperties operator+(const RemovalProperties& b) const
    //    {
    //        return RemovalProperties{mass + b.mass,
    //                                 {
    //                                     momentum[0] + b.momentum[0],
    //                                     momentum[1] + b.momentum[1],
    //                                     momentum[2] + b.momentum[2],
    //                                 },
    //                                 count + b.count};
    //    }
};

__device__ void atomicAdd(RemovalProperties* x, const RemovalProperties& y)
{
    atomicAdd(&(x->mass), y.mass);
    atomicAdd(&(x->momentum[0]), y.momentum[0]);
    atomicAdd(&(x->momentum[1]), y.momentum[1]);
    atomicAdd(&(x->momentum[2]), y.momentum[2]);
    atomicAdd(&(x->count), y.count);
}

static __device__ RemovalProperties device_accreted;
static __device__ RemovalProperties device_removed;

using cstone::TravConfig;

template<typename T1, typename Th, typename Tkeys, typename T2, typename Tm, typename Tv>
__global__ void computeAccretionConditionKernel(size_t first, size_t last, const T1* x, const T1* y, const T1* z,
                                                const Th* h, Tkeys* keys, const Tm* m, const Tv* vx, const Tv* vy,
                                                const Tv* vz, T2 star_x, T2 star_y, T2 star_z, T2 star_size2,
                                                T2 removal_limit_h)
{
    cstone::LocalIndex i = first + blockDim.x * blockIdx.x + threadIdx.x;

    RemovalProperties accreted{};
    RemovalProperties removed{};

    if (i >= last) {}
    else
    {
        const double dx    = x[i] - star_x;
        const double dy    = y[i] - star_y;
        const double dz    = z[i] - star_z;
        const double dist2 = dx * dx + dy * dy + dz * dz;

        if (dist2 < star_size2)
        {
            // Accrete on star
            keys[i]              = cstone::removeKey<Tkeys>::value;
            accreted.mass        = m[i];
            accreted.momentum[0] = m[i] * vx[i];
            accreted.momentum[1] = m[i] * vy[i];
            accreted.momentum[2] = m[i] * vz[i];
            accreted.count       = 1;
        }
        else if (h[i] > removal_limit_h)
        {
            // Remove from system
            keys[i]             = cstone::removeKey<Tkeys>::value;
            removed.mass        = m[i];
            removed.momentum[0] = m[i] * vx[i];
            removed.momentum[1] = m[i] * vy[i];
            removed.momentum[2] = m[i] * vz[i];
            removed.count       = 1;
        }
    }

    typedef cub::BlockReduce<RemovalProperties, TravConfig::numThreads> BlockReduce;
    __shared__ typename BlockReduce::TempStorage                        temp_accreted, temp_removed;

    struct Sum
    {
        __device__ RemovalProperties operator()(const RemovalProperties& a, const RemovalProperties& b) const
        {
            return RemovalProperties{a.mass + b.mass,
                                     {
                                         a.momentum[0] + b.momentum[0],
                                         a.momentum[1] + b.momentum[1],
                                         a.momentum[2] + b.momentum[2],
                                     },
                                     a.count + b.count};
        }
    };

    RemovalProperties block_accreted = BlockReduce(temp_accreted).Reduce(accreted, Sum{});
    RemovalProperties block_removed  = BlockReduce(temp_removed).Reduce(removed, Sum{});

    __syncthreads();

    if (threadIdx.x == 0)
    {
        atomicAdd(&device_accreted, block_accreted);
        atomicAdd(&device_removed, block_removed);
    }
}

template<typename Dataset, typename StarData>
void computeAccretionConditionGPU(size_t first, size_t last, Dataset& d, StarData& star)
{
    cstone::LocalIndex numParticles = last - first;
    unsigned           numThreads   = 256;
    unsigned           numBlocks    = (numParticles + numThreads - 1) / numThreads;

    RemovalProperties accreted_local{}, removed_local{};

    cudaMemcpyToSymbol(device_accreted, &accreted_local, sizeof(accreted_local));
    cudaMemcpyToSymbol(device_removed, &removed_local, sizeof(removed_local));

    computeAccretionConditionKernel<<<numBlocks, numThreads>>>(
        first, last, rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z), rawPtr(d.devData.h),
        rawPtr(d.devData.keys), rawPtr(d.devData.m), rawPtr(d.devData.vx), rawPtr(d.devData.vy), rawPtr(d.devData.vz),
        star.position[0], star.position[1], star.position[2], star.inner_size * star.inner_size, star.removal_limit_h);
    checkGpuErrors(cudaGetLastError());
    checkGpuErrors(cudaDeviceSynchronize());

    cudaMemcpyFromSymbol(&accreted_local, device_accreted, sizeof(accreted_local));
    cudaMemcpyFromSymbol(&removed_local, device_removed, sizeof(removed_local));

    star.m_accreted_local    = accreted_local.mass;
    star.p_accreted_local[0] = accreted_local.momentum[0];
    star.p_accreted_local[1] = accreted_local.momentum[1];
    star.p_accreted_local[2] = accreted_local.momentum[2];
    star.n_accreted_local    = accreted_local.count;

    star.m_removed_local    = removed_local.mass;
    star.p_removed_local[0] = removed_local.momentum[0];
    star.p_removed_local[1] = removed_local.momentum[1];
    star.p_removed_local[2] = removed_local.momentum[2];
    star.n_removed_local    = removed_local.count;
}

template void computeAccretionConditionGPU(size_t, size_t, sphexa::ParticlesData<cstone::GpuTag>&, StarData&);
