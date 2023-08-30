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

    auto dataTuple()
    {
        auto ret = std::tie(hydro, chem);
        return ret;
    }

    auto data()
    {
        using FieldVariant = std::variant<HydroData*, ChemData*>;
        return std::apply([](auto&... fields) { return std::array<FieldVariant, sizeof...(fields)>{&fields...}; },
                          dataTuple());
    }

    void setOutputFields(const std::vector<std::string>& outFields)
    {
        const auto&                                      ds         = data();
        constexpr size_t                                 n_datasets = std::tuple_size_v<decltype(dataTuple())>;
        std::array<std::vector<std::string>, n_datasets> outFieldsFwd;

        auto assign_if_correct =
            [](const auto* d, const std::string& outFieldName, std::vector<std::string>& outFieldsVec)
        {
            std::string prefix(d->datasetPrefix);
            const auto  m = std::mismatch(prefix.begin(), prefix.end(), outFieldName.begin());
            if (m.first == prefix.end())
            {
                outFieldsVec.emplace_back(m.second, outFieldName.end());
                return true;
            }
            return false;
        };

        auto assign_field_to_dataset = [&](const std::string& outFieldName)
        {
            for (size_t i = outFieldsFwd.size() - 1; i >= 0; i--)
            {
                const bool correct = std::visit(
                    [&](const auto* d) { return assign_if_correct(d, outFieldName, outFieldsFwd[i]); }, ds[i]);
                if (correct) { return; }
            }
        };
        std::for_each(outFields.begin(), outFields.end(), assign_field_to_dataset);

        for (size_t i = 0; i < n_datasets; i++)
        {
            std::visit([&](auto* d) { d->setOutputFields(outFieldsFwd[i]); }, ds[i]);
        }
    }
};

} // namespace sphexa
