#pragma once

#include "cstone/findneighbors.hpp"

namespace sph
{

using cstone::LocalIndex;

template<class T, class KeyType>
void findNeighborsSph(const T* x, const T* y, const T* z, T* h, LocalIndex firstId, LocalIndex lastId,
                      const cstone::Box<T>& box, const cstone::OctreeNsView<T, KeyType>& treeView, unsigned ng0,
                      unsigned ngmax, LocalIndex* neighbors, unsigned* nc)
{
    LocalIndex numWork = lastId - firstId;

    // unsigned ngmin = ng0 / 4;
    unsigned ngmin = ng0;

    size_t        numFails     = 0;
    constexpr int maxIteration = 1000;

#pragma omp parallel for reduction(+ : numFails)
    for (LocalIndex i = 0; i < numWork; ++i)
    {
        LocalIndex id    = i + firstId;
        unsigned   ncSph = 1 + findNeighbors(id, x, y, z, h, treeView, box, ngmax, neighbors + i * ngmax);

        int iteration = 0;
        while (ngmin > ncSph || (ncSph - 1) > ngmax && iteration++ < maxIteration)
        {
            h[id] = updateH(ng0, ncSph, h[id]);
            ncSph = 1 + findNeighbors(id, x, y, z, h, treeView, box, ngmax, neighbors + i * ngmax);
        }
        numFails += (iteration == maxIteration);

        nc[i] = ncSph;
    }

    // exact neighbour
#pragma omp parallel for
    for (LocalIndex i = 0; i < numWork; ++i)
    {
        unsigned n = nc[i];
        if (n < ng0) continue;
        std::vector<T>      dist(n);
        std::vector<size_t> arg(n);
        for (size_t j = 0; j < n; j++)
        {
            unsigned nb_j = *(neighbors + i * ngmax + j);
            dist[j]       = (x[nb_j] * x[nb_j] + y[nb_j] * y[nb_j] + z[nb_j] * z[nb_j]);
        }
        std::iota(arg.begin(), arg.end(), 0l);
        std::sort(arg.begin(), arg.end(), [&dist](size_t i, size_t j) { return dist[i] > dist[j]; });

        std::vector<size_t> sorted(n);
        for (size_t j = 0; j < n; j++)
        {
            sorted[j] = *(neighbors + i * ngmax + arg[j]);
        }
        for (size_t j = 0; j < n; j++)
        {
            *(neighbors + i * ngmax + j) = sorted[j];
        }
    }
        if (numFails)
        {
            std::cout << "Coupled h-neighbor count updated failed to converge for " << numFails << " particles"
                      << std::endl;
        }
    }

    //! @brief perform neighbor search together with updating the smoothing lengths
    template<class T, class Dataset>
    void findNeighborsSfc(size_t startIndex, size_t endIndex, Dataset & d, const cstone::Box<T>& box)
    {
        if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{}) { return; }

        if (d.ng0 > d.ngmax) { throw std::runtime_error("ng0 should be smaller than ngmax\n"); }

        findNeighborsSph(d.x.data(), d.y.data(), d.z.data(), d.h.data(), startIndex, endIndex, box, d.treeView.nsView(),
                         d.ng0, d.ngmax, d.neighbors.data(), d.nc.data() + startIndex);
    }

} // namespace sph
