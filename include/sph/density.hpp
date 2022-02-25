#pragma once

#include <vector>

#include "cstone/findneighbors.hpp"

#include "kernels.hpp"
#include "task.hpp"
#include "kernel/computeDensity.hpp"
#ifdef USE_CUDA
#include "cuda/sph.cuh"
#endif

namespace sphexa
{
namespace sph
{
template<typename T, class Dataset>
void computeDensityImpl(const Task& t, Dataset& d, const cstone::Box<T>& box)
{
    // number of particles in task
    size_t numParticles = t.size();

    const size_t ngmax          = t.ngmax;
    const int*   neighbors      = t.neighbors.data();
    const int*   neighborsCount = t.neighborsCount.data();

    const T* h = d.h.data();
    const T* m = d.m.data();
    const T* x = d.x.data();
    const T* y = d.y.data();
    const T* z = d.z.data();

    const T* wh  = d.wh.data();
    const T* whd = d.whd.data();

    T* ro = d.ro.data();

    const T K         = d.K;
    const T sincIndex = d.sincIndex;

#if defined(USE_OMP_TARGET)
    const size_t np           = d.x.size();
    const size_t ltsize       = d.wh.size();
    const size_t n            = numParticles;
    const size_t allNeighbors = n * ngmax;

// clang-format off
#pragma omp target map(to                                                                                                                  \
                       : n, neighbors [:allNeighbors], neighborsCount [:n], m [0:np], h [0:np], x [0:np], y [0:np], z [0:np],  wh [0:ltsize], whd [0:ltsize])    \
                   map(from                                                                                                                \
                       : ro [:n])
#pragma omp teams distribute parallel for
// clang-format on
#else
#pragma omp parallel for
#endif
    for (size_t pi = 0; pi < numParticles; pi++)
    {
        // int neighLoc[ngmax];
        // int count;
        // cstone::findNeighbors(
        //    pi, x, y, z, h, box, cstone::sfcKindPointer(d.codes.data()), neighLoc, &count, d.codes.size(), ngmax);

        int i = pi + t.firstParticle;

        ro[i] = kernels::densityJLoop(
            i, sincIndex, K, box, neighbors + ngmax * pi, neighborsCount[pi], x, y, z, h, m, wh, whd);

#ifndef NDEBUG
        if (std::isnan(ro[i]))
            printf("ERROR::Density(%zu) density %f, position: (%f %f %f), h: %f\n", pi, ro[i], x[i], y[i], z[i], h[i]);
#endif
    }
}

template<typename T, class Dataset>
void computeDensity(std::vector<Task>& taskList, Dataset& d, const cstone::Box<T>& box)
{
#if defined(USE_CUDA)
    cuda::computeDensity<Dataset>(taskList, d, box);
#else
    for (const auto& task : taskList)
    {
        computeDensityImpl<T>(task, d, box);
    }
#endif
}

} // namespace sph
} // namespace sphexa
