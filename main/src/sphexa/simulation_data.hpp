/*
 * MIT License
 *
 * Copyright (c) 2022 CSCS, ETH Zurich, University of Basel, University of Zurich
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
 * @brief Contains the object holding all simulation data
 *
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

#pragma once

#include <mpi.h>

#include "cooling/chemistry_data.hpp"
#include "sph/particles_data.hpp"

namespace sphexa
{

//! @brief the place to store hydro, chemistry, nuclear and other simulation data
template<typename T, typename KeyType_, class AccType>
class SimulationData
{
public:
    using AcceleratorType = AccType;
    using KeyType         = KeyType_;
    using RealType        = T;

    using HydroData = ParticlesData<RealType, KeyType, AccType>;
    using ChemData  = cooling::ChemistryData<T>;

    //! @brief spacially distributed data for hydrodynamics and gravity
    HydroData hydro;

    //! @brief chemistry data for radiative cooling, e.g. for GRACKLE
    ChemData chem;

    //! @brief non-spacially distributed nuclear abundances
    // NuclearData nuclear;

    MPI_Comm comm;

    inline static constexpr std::array dataPrefix{"", "chem::"};

    auto dataTuple()
    {
        auto ret = std::tie(hydro, chem);
        static_assert(std::tuple_size_v<decltype(ret)> == dataPrefix.size());
        return ret;
    }

    void setOutputFields(const std::vector<std::string>& outFields)
    {
        std::array<std::vector<std::string>, dataPrefix.size()> outFieldsFwd;

        auto assign_name = [&outFieldsFwd](std::string outField)
        {
            for (size_t i = dataPrefix.size() - 1; i >= 0; i--)
            {
                std::string prefix(dataPrefix[i]);
                const auto m = std::mismatch(prefix.begin(), prefix.end(), outField.begin());
                if (m.first == prefix.end())
                {
                    outFieldsFwd[i].emplace_back(m.second, outField.end());
                    break;
                }
            }
        };

        std::for_each(outFields.begin(), outFields.end(), assign_name);

        auto data              = dataTuple();
        auto outFieldsFwdTuple = std::tuple_cat(outFieldsFwd);
        for_each_tuple_zip([](auto& d, const std::vector<std::string>& f) { d.setOutputFields(f); }, data,
                           outFieldsFwdTuple);
    }
};

} // namespace sphexa
