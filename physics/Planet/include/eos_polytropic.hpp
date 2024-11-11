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
 * @brief Density i-loop OpenMP driver
 *
 * @author Ruben Cabezon <ruben.cabezon@unibas.ch>
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

#pragma once

#include "cstone/cuda/cuda_utils.hpp"
#include "eos_polytropic_gpu.hpp"
#include "eos_polytropic_loop.hpp"
#include "sph/particles_data_stubs.hpp"
#include "sph/eos.hpp"

namespace planet
{

template<typename Dataset, typename StarData>
void computePolytropic_HydroStdImpl(size_t startIndex, size_t endIndex, Dataset& d, const StarData& star)
{
    const auto* rho = d.rho.data();

    auto* p = d.p.data();
    auto* c = d.c.data();

#pragma omp parallel for schedule(static)
    for (size_t i = startIndex; i < endIndex; ++i)
    {
        std::tie(p[i], c[i]) = polytropicEOS(star.Kpoly, star.exp_poly, d.gamma, rho[i]);
    }
}

template<class Dataset, typename StarData>
void computePolytropicEOS_HydroStd(size_t startIndex, size_t endIndex, Dataset& d, const StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        computePolytropicEOS_HydroStdGPU(startIndex, endIndex, star.Kpoly, star.exp_poly, d.gamma,
                                            rawPtr(d.devData.rho), rawPtr(d.devData.p), rawPtr(d.devData.c));
    }
    else { computePolytropic_HydroStdImpl(startIndex, endIndex, d, star); }
}

} // namespace sph
