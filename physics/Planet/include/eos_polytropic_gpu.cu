/*
 * MIT License
 *
 * Copyright (c) 2021 CSCS, ETH Zurich
 *               2021 University of Basel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*! @file
 * @brief Density i-loop GPU driver
 *
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

//#include "cstone/cuda/cuda_utils.cuh"
#include "cstone/primitives/math.hpp"
#include "cstone/util/tuple.hpp"
#include "eos_polytropic_gpu.hpp"
#include "eos_polytropic_loop.hpp"
#include "sph/particles_data.hpp"
#include "star_data.hpp"

namespace planet
{

template<typename T1, typename T2, typename T3, typename Trho, typename Tp, typename Tc>
__global__ void computePolytropicEOS_HydroStdKernel(size_t firstParticle, size_t lastParticle, T1 Kpoly, T2 exp_poly,
                                                    T3 gamma, const Trho* rho, Tp* p, Tc* c)
{
    unsigned i = firstParticle + blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= lastParticle) return;

    util::tie(p[i], c[i]) = polytropicEOS(Kpoly, exp_poly, gamma, rho[i]);
}

template<typename Dataset, typename StarData>
void computePolytropicEOS_HydroStdGPU(size_t firstParticle, size_t lastParticle, Dataset& d, const StarData& star)
{
    if (firstParticle == lastParticle) { return; }
    unsigned numThreads = 256;
    unsigned numBlocks  = cstone::iceil(lastParticle - firstParticle, numThreads);
    computePolytropicEOS_HydroStdKernel<<<numBlocks, numThreads>>>(firstParticle, lastParticle, star.Kpoly,
                                                                   star.exp_poly, d.gamma, rawPtr(d.devData.rho),
                                                                   rawPtr(d.devData.p), rawPtr(d.devData.c));

    checkGpuErrors(cudaDeviceSynchronize());
}

template void computePolytropicEOS_HydroStdGPU(size_t, size_t, sphexa::ParticlesData<cstone::GpuTag>&, const StarData&);
} // namespace planet
