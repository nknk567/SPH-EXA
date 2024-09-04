//
// Created by Noah Kubli on 28.08.2024.
//

#pragma once

#include <algorithm>

#include "cstone/cuda/gpu_config.cuh"
#include "cstone/cuda/cuda_utils.cuh"
#include "cstone/primitives/warpscan.cuh"
#include "cstone/sfc/box.hpp"
#include "cstone/findneighbors.hpp"

#include "cstone/traversal/find_neighbors.cuh"

#include "cstone/tree/definitions.h"
#include "cstone/tree/octree.hpp"

namespace cstone
{

/*! @brief count neighbors within a cutoff
 *
 * @tparam       T             float or double
 * @param[in]    sourceBody    source body x,y,z
 * @param[in]    validLaneMask number of lanes that contain valid source bodies
 * @param[in]    pos_i         target body x,y,z,h
 * @param[in]    box           global coordinate bounding box
 * @param[in]    sourceBodyIdx index of source body of each lane
 * @param[in]    ngmax         maximum number of neighbors to store
 * @param[inout] nc_i          target neighbor count to add to
 * @param[out]   nidx_i        indices of neighbors of the target body to store
 *
 * Number of computed particle-particle pairs per call is GpuConfig::warpSize^2 * TravConfig::nwt
 */

__device__ unsigned left_child(unsigned i) { return 2 * i + 1; }
__device__ unsigned right_child(unsigned i) { return 2 * i + 2; }
__device__ unsigned parent(unsigned i) { return (i - 1) / 2; }

template<typename T>
__device__ void swap_values(T& a, T& b)
{
    T temp{a};
    a = b;
    b = temp;
}

template<unsigned stride>
__device__ void push(unsigned* nidx, double* d2, unsigned nc_i)
{
    if (nc_i <= 1) return;

    unsigned i = nc_i - 1;
    while (i > 0 && d2[stride * parent(i)] < d2[stride * i])
    {
        swap_values(nidx[stride * parent(i)], nidx[stride * i]);
        swap_values(d2[stride * parent(i)], d2[stride * i]);
        i = parent(i);
    }
}

template<unsigned stride>
__device__ void pop(unsigned* nidx, double* d2, unsigned nc_i)
{
    if (nc_i <= 1) return;

    swap_values(nidx[stride * (nc_i - 1)], nidx[stride * 0]);
    swap_values(d2[stride * (nc_i - 1)], d2[stride * 0]);
    nc_i--;

    unsigned i     = 0;
    unsigned left  = left_child(i);
    unsigned right = right_child(i);
    while (true)
    {
        unsigned new_index = i;
        if (left < nc_i && d2[stride * left] > d2[stride * i]) { new_index = left; }
        if (right < nc_i && d2[stride * right] > d2[stride * new_index]) { new_index = right; }
        // if (i == new_index) break; //balloc-sync
        if (!cstone::ballotSync(i != new_index)) { break; }
        swap_values(nidx[stride * new_index], nidx[stride * i]);
        swap_values(d2[stride * new_index], d2[stride * i]);

        left  = left_child(new_index);
        right = right_child(new_index);
        i     = new_index;
    }
}

template<bool UsePbc, class Tc>
__device__ void countNeighborsHeap(Vec3<Tc> sourceBody,
                                   int numLanesValid,
                                   const util::array<Vec4<Tc>, TravConfig::nwt>& pos_i,
                                   const Box<Tc>& box,
                                   cstone::LocalIndex sourceBodyIdx,
                                   unsigned ngmax,
                                   unsigned nc_i[TravConfig::nwt],
                                   unsigned* nidx_i,
                                   double* d2_i)
{
    unsigned laneIdx = threadIdx.x & (GpuConfig::warpSize - 1);

    for (int j = 0; j < numLanesValid; j++)
    {
        Vec3<Tc> pos_j{shflSync(sourceBody[0], j), shflSync(sourceBody[1], j), shflSync(sourceBody[2], j)};
        cstone::LocalIndex idx_j = shflSync(sourceBodyIdx, j);

#pragma unroll
        for (int k = 0; k < TravConfig::nwt; k++)
        {
            Tc d2 = distanceSq<UsePbc>(pos_j[0], pos_j[1], pos_j[2], pos_i[k][0], pos_i[k][1], pos_i[k][2], box);
            if (d2 < pos_i[k][3] && d2 > Tc(0.0))
            {
                unsigned first_k = laneIdx + k * GpuConfig::warpSize;
                if (nc_i[k] < ngmax)
                {
                    nidx_i[nc_i[k] * TravConfig::targetSize + first_k] = idx_j;
                    d2_i[nc_i[k] * TravConfig::targetSize + first_k]   = d2;

                    push<TravConfig::targetSize>(nidx_i + first_k, d2_i + first_k, nc_i[k]);
                    nc_i[k]++;
                }
                else if (d2 < d2_i[(0) * TravConfig::targetSize + first_k])
                {
                    pop<TravConfig::targetSize>(nidx_i + first_k, d2_i + first_k, nc_i[k]);

                    nidx_i[(nc_i[k] - 1) * TravConfig::targetSize + first_k] = idx_j;
                    d2_i[(nc_i[k] - 1) * TravConfig::targetSize + first_k]   = d2;

                    push<TravConfig::targetSize>(nidx_i + first_k, d2_i + first_k, nc_i[k]);
                }
            }
        }
    }
}

/*! @brief traverse one warp with up to TravConfig::targetSize target bodies down the tree
 *
 * @param[inout] nc_i           output neighbor counts to add to, TravConfig::nwt per lane
 * @param[in]    nidx_i         output neighbor indices, up to @p ngmax * TravConfig::nwt per lane will be stored
 * @param[in]    ngmax          maximum number of neighbor particle indices to store in @p nidx_i
 * @param[in]    pos_i          target x,y,z,4h^2, TravConfig::nwt per lane
 * @param[in]    targetCenter   geometrical target center
 * @param[in]    targetSize     geometrical target bounding box size
 * @param[in]    x,y,z,h        source bodies as referenced by tree cells
 * @param[in]    tree           octree data view
 * @param[in]    initNodeIdx    traversal will be started with all children of the parent of @p initNodeIdx
 * @param[in]    box            global coordinate bounding box
 * @param[-]     tempQueue      shared mem int pointer to GpuConfig::warpSize ints, uninitialized
 * @param[-]     cellQueue      pointer to global memory, size defined by TravConfig::memPerWarp, uninitialized
 * @return                      Number of P2P interactions tested to the group of target particles.
 *                              The total for the warp is the numbers returned here times the number of valid
 *                              targets in the warp.
 *
 * Constant input pointers are additionally marked __restrict__ to indicate to the compiler that loads
 * can be routed through the read-only/texture cache.
 */
template<bool UsePbc, class Tc, class Th, class KeyType>
__device__ uint2 traverseWarpHeap(unsigned* nc_i,
                                  unsigned* nidx_i,
                                  double* d2_i,
                                  unsigned ngmax,
                                  const util::array<Vec4<Tc>, TravConfig::nwt>& pos_i,
                                  const Vec3<Tc> targetCenter,
                                  const Vec3<Tc> targetSize,
                                  const Tc* __restrict__ x,
                                  const Tc* __restrict__ y,
                                  const Tc* __restrict__ z,
                                  const Th* __restrict__ /*h*/,
                                  const OctreeNsView<Tc, KeyType>& tree,
                                  int initNodeIdx,
                                  const Box<Tc>& box,
                                  volatile int* tempQueue,
                                  int* cellQueue)
{
    const TreeNodeIndex* __restrict__ childOffsets   = tree.childOffsets;
    const TreeNodeIndex* __restrict__ internalToLeaf = tree.internalToLeaf;
    const LocalIndex* __restrict__ layout            = tree.layout;
    const Vec3<Tc>* __restrict__ centers             = tree.centers;
    const Vec3<Tc>* __restrict__ sizes               = tree.sizes;

    const int laneIdx = threadIdx.x & (GpuConfig::warpSize - 1);

    unsigned p2pCounter = 0, maxStack = 0;

    int bodyQueue; // warp queue for source body indices

    // populate initial cell queue
    if (laneIdx == 0) { cellQueue[0] = initNodeIdx; }

    // these variables are always identical on all warp lanes
    int numSources   = 1; // current stack size
    int newSources   = 0; // stack size for next level
    int oldSources   = 0; // cell indices done
    int sourceOffset = 0; // current level stack pointer, once this reaches numSources, the level is done
    int bdyFillLevel = 0; // fill level of the source body warp queue

    while (numSources > 0) // While there are source cells to traverse
    {
        int sourceIdx   = sourceOffset + laneIdx; // Source cell index of current lane
        int sourceQueue = 0;
        if (laneIdx < GpuConfig::warpSize / 8)
        {
            sourceQueue = cellQueue[ringAddr(oldSources + sourceIdx)]; // Global source cell index in queue
        }
        sourceQueue = spreadSeg8(sourceQueue);
        sourceIdx   = shflSync(sourceIdx, laneIdx >> 3);

        const Vec3<Tc> curSrcCenter = centers[sourceQueue];      // Current source cell center
        const Vec3<Tc> curSrcSize   = sizes[sourceQueue];        // Current source cell center
        const int childBegin        = childOffsets[sourceQueue]; // First child cell
        const bool isNode           = childBegin;
        const bool isClose          = cellOverlap<UsePbc>(curSrcCenter, curSrcSize, targetCenter, targetSize, box);
        const bool isSource         = sourceIdx < numSources; // Source index is within bounds
        const bool isDirect         = isClose && !isNode && isSource;
        const int leafIdx           = (isDirect) ? internalToLeaf[sourceQueue] : 0; // the cstone leaf index

        // Split
        const bool isSplit     = isNode && isClose && isSource;                   // Source cell must be split
        const int numChildLane = exclusiveScanBool(isSplit);                      // Exclusive scan of numChild
        const int numChildWarp = reduceBool(isSplit);                             // Total numChild of current warp
        sourceOffset += imin(GpuConfig::warpSize / 8, numSources - sourceOffset); // advance current level stack pointer
        int childIdx = oldSources + numSources + newSources + numChildLane;       // Child index of current lane
        if (isSplit) { cellQueue[ringAddr(childIdx)] = childBegin; }              // Queue child cells for next level
        newSources += numChildWarp; // Increment source cell count for next loop

        // check for cellQueue overflow
        const unsigned stackUsed = newSources + numSources - sourceOffset; // current cellQueue size
        maxStack                 = max(stackUsed, maxStack);
        if (stackUsed > TravConfig::memPerWarp) { return {0xFFFFFFFF, maxStack}; } // Exit if cellQueue overflows

        // Direct
        const int firstBody     = layout[leafIdx];
        const int numBodies     = (layout[leafIdx + 1] - firstBody) & -int(isDirect); // Number of bodies in cell
        bool directTodo         = numBodies;
        const int numBodiesScan = inclusiveScanInt(numBodies);                      // Inclusive scan of numBodies
        int numBodiesLane       = numBodiesScan - numBodies;                        // Exclusive scan of numBodies
        int numBodiesWarp       = shflSync(numBodiesScan, GpuConfig::warpSize - 1); // Total numBodies of current warp
        int prevBodyIdx         = 0;
        while (numBodiesWarp > 0) // While there are bodies to process from current source cell set
        {
            tempQueue[laneIdx] = 1; // Default scan input is 1, such that consecutive lanes load consecutive bodies
            if (directTodo && (numBodiesLane < GpuConfig::warpSize))
            {
                directTodo               = false;          // Set cell as processed
                tempQueue[numBodiesLane] = -1 - firstBody; // Put first source cell body index into the queue
            }
            const int bodyIdx = inclusiveSegscanInt(tempQueue[laneIdx], prevBodyIdx);
            // broadcast last processed bodyIdx from the last lane to restart the scan in the next iteration
            prevBodyIdx = shflSync(bodyIdx, GpuConfig::warpSize - 1);

            if (numBodiesWarp >= GpuConfig::warpSize) // Process bodies from current set of source cells
            {
                // Load source body coordinates
                const Vec3<Tc> sourceBody = {x[bodyIdx], y[bodyIdx], z[bodyIdx]};
                countNeighborsHeap<UsePbc>(sourceBody, GpuConfig::warpSize, pos_i, box, bodyIdx, ngmax, nc_i, nidx_i,
                                           d2_i);
                numBodiesWarp -= GpuConfig::warpSize;
                numBodiesLane -= GpuConfig::warpSize;
                p2pCounter += GpuConfig::warpSize;
            }
            else // Fewer than warpSize bodies remaining from current source cell set
            {
                // push the remaining bodies into bodyQueue
                int topUp = shflUpSync(bodyIdx, bdyFillLevel);
                bodyQueue = (laneIdx < bdyFillLevel) ? bodyQueue : topUp;

                bdyFillLevel += numBodiesWarp;
                if (bdyFillLevel >= GpuConfig::warpSize) // If this causes bodyQueue to spill
                {
                    // Load source body coordinates
                    const Vec3<Tc> sourceBody = {x[bodyQueue], y[bodyQueue], z[bodyQueue]};
                    countNeighborsHeap<UsePbc>(sourceBody, GpuConfig::warpSize, pos_i, box, bodyQueue, ngmax, nc_i,
                                               nidx_i, d2_i);
                    bdyFillLevel -= GpuConfig::warpSize;
                    // bodyQueue is now empty; put body indices that spilled into the queue
                    bodyQueue = shflDownSync(bodyIdx, numBodiesWarp - bdyFillLevel);
                    p2pCounter += GpuConfig::warpSize;
                }
                numBodiesWarp = 0; // No more bodies to process from current source cells
            }
        }

        //  If the current level is done
        if (sourceOffset >= numSources)
        {
            oldSources += numSources;      // Update finished source size
            numSources   = newSources;     // Update current source size
            sourceOffset = newSources = 0; // Initialize next source size and offset
        }
    }

    if (bdyFillLevel > 0) // If there are leftover direct bodies
    {
        const bool laneHasBody = laneIdx < bdyFillLevel;
        // Load position of source bodies, with padding for invalid lanes
        const Vec3<Tc> sourceBody =
            laneHasBody ? Vec3<Tc>{x[bodyQueue], y[bodyQueue], z[bodyQueue]} : Vec3<Tc>{Tc(0), Tc(0), Tc(0)};
        countNeighborsHeap<UsePbc>(sourceBody, bdyFillLevel, pos_i, box, bodyQueue, ngmax, nc_i, nidx_i, d2_i);
        p2pCounter += bdyFillLevel;
    }

    return {p2pCounter, maxStack};
}

/*! @brief Find neighbors of a group of given particles, does not count self reference: min return value is 0
 *
 * @param[in]  bodyBegin   index of first particle in (x,y,z) to look for neighbors
 * @param[in]  bodyEnd     last (excluding) index of particle to look for neighbors
 * @param[in]  x           particle x coordinates
 * @param[in]  y           particle y coordinates
 * @param[in]  z           particle z coordinates
 * @param[in]  h           particle smoothing lengths
 * @param[in]  tree        octree connectivity and cell data
 * @param[in]  box         global coordinate bounding box
 * @param[out] warpNidx    storage for up to ngmax neighbor part. indices for each of the (bodyEnd - bodyBegin) targets
 * @param[in]  ngmax       maximum number of neighbors to store
 * @param[-]   globalPool  global memory for cell traversal stack
 * @return                 actual neighbor count of the particle handled by the executing warp lane, can be > ngmax,
 *                         minimum returned value is 0
 *
 * Note: Number of handled particles (bodyEnd - bodyBegin) should be GpuConfig::warpSize * TravConfig::nwt or smaller
 */
template<class Tc, class Th, class KeyType>
__device__ util::array<unsigned, TravConfig::nwt> traverseNeighborsHeap(cstone::LocalIndex bodyBegin,
                                                                        cstone::LocalIndex bodyEnd,
                                                                        const Tc* __restrict__ x,
                                                                        const Tc* __restrict__ y,
                                                                        const Tc* __restrict__ z,
                                                                        const Th* __restrict__ h,
                                                                        const OctreeNsView<Tc, KeyType>& tree,
                                                                        const Box<Tc>& box,
                                                                        cstone::LocalIndex* warpNidx,
                                                                        double* warpd2,
                                                                        unsigned ngmax,
                                                                        TreeNodeIndex* globalPool)
{
    const unsigned laneIdx = threadIdx.x & (GpuConfig::warpSize - 1);
    const unsigned warpIdx = threadIdx.x >> GpuConfig::warpSizeLog2;

    constexpr unsigned numWarpsPerBlock = TravConfig::numThreads / GpuConfig::warpSize;

    __shared__ int sharedPool[TravConfig::numThreads];

    // warp-common shared mem, 1 int per thread
    int* tempQueue = sharedPool + GpuConfig::warpSize * warpIdx;
    // warp-common global mem storage
    int* cellQueue = globalPool + TravConfig::memPerWarp * ((blockIdx.x * numWarpsPerBlock) + warpIdx);

    util::array<Vec4<Tc>, TravConfig::nwt> pos_i = loadTarget(bodyBegin, bodyEnd, laneIdx, x, y, z, h);
    auto [targetCenter, targetSize]              = warpBbox(pos_i);
    targetSize *= Tc(tree.searchExtFactor);

#pragma unroll
    for (int k = 0; k < TravConfig::nwt; ++k)
    {
        auto r      = pos_i[k][3];
        pos_i[k][3] = r * r;
    }

    auto pbc    = BoundaryType::periodic;
    bool anyPbc = box.boundaryX() == pbc || box.boundaryY() == pbc || box.boundaryZ() == pbc;
    bool usePbc = anyPbc && !insideBox(targetCenter, targetSize, box);

    util::array<unsigned, TravConfig::nwt> nc_i; // NOLINT
    nc_i = 0;

    // start traversal with node 1 (first child of the root), implies siblings as well
    // if traversal should be started at node x, then initNode should be set to the first child of x
    int initNode = 1;

    uint2 warpStats;
    if (usePbc)
    {
        warpStats = traverseWarpHeap<true>(nc_i.data(), warpNidx, warpd2, ngmax, pos_i, targetCenter, targetSize, x, y,
                                           z, h, tree, initNode, box, tempQueue, cellQueue);
    }
    else
    {
        warpStats = traverseWarpHeap<false>(nc_i.data(), warpNidx, warpd2, ngmax, pos_i, targetCenter, targetSize, x, y,
                                            z, h, tree, initNode, box, tempQueue, cellQueue);
    }
    unsigned numP2P   = warpStats.x;
    unsigned maxStack = warpStats.y;
    assert(numP2P != 0xFFFFFFFF);

    if (laneIdx == 0)
    {
        unsigned targetGroupSize = bodyEnd - bodyBegin;
        atomicAdd(&ncStats[NcStats::sumP2P], NcStats::type(numP2P) * targetGroupSize);
        atomicMax(&ncStats[NcStats::maxP2P], NcStats::type(numP2P));
        atomicMax(&ncStats[NcStats::maxStack], NcStats::type(maxStack));
    }

    return nc_i;
}

} // namespace cstone
