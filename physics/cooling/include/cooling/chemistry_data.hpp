/*
 * MIT License
 *
 * Copyright (c) 2022 CSCS, ETH Zurich, University of Zurich, University of Basel
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
 * @brief Unit tests for ParticlesData
 *
 * @author Noah Kubli <noah.kubli@uzh.ch>
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */
#pragma once

#include <array>
#include <optional>

#include "cstone/fields/enumerate.hpp"
#include "cstone/fields/field_states.hpp"
#include "cstone/fields/particles_get.hpp"
#include "cstone/util/reallocate.hpp"
#include "cstone/util/traits.hpp"
#include "cooling/cooler.hpp"

namespace cooling
{

template<class T>
class ChemistryData : public cstone::FieldStates<ChemistryData<T>>
{
public:
    inline static constexpr size_t numFields = 21;

    template<class ValueType>
    using FieldVector     = std::vector<ValueType, std::allocator<ValueType>>;
    using RealType        = T;
    using AcceleratorType = cstone::CpuTag;
    using FieldVariant =
        std::variant<FieldVector<float>*, FieldVector<double>*, FieldVector<unsigned>*, FieldVector<uint64_t>*>;

    std::array<FieldVector<T>, numFields> fields;

    auto dataTuple() { return dataTuple_helper(std::make_index_sequence<numFields>{}); }

    auto data()
    {
        return std::apply([](auto&... fields) { return std::array<FieldVariant, sizeof...(fields)>{&fields...}; },
                          dataTuple());
    }

    void resize(size_t size)
    {
        double growthRate = 1.05;
        auto   data_      = data();

        for (size_t i = 0; i < data_.size(); ++i)
        {
            if (this->isAllocated(i))
            {
                std::visit([size, growthRate](auto& arg) { reallocate(*arg, size, growthRate); }, data_[i]);
            }
        }
    }

    //! Grackle field names
    inline static constexpr std::array fieldNames{"HI_fraction",
                                                  "HII_fraction",
                                                  "HM_fraction",
                                                  "HeI_fraction",
                                                  "HeII_fraction",
                                                  "HeIII_fraction",
                                                  "H2I_fraction",
                                                  "H2II_fraction",
                                                  "DI_fraction",
                                                  "DII_fraction",
                                                  "HDI_fraction",
                                                  "e_fraction",
                                                  "metal_fraction",
                                                  "volumetric_heating_rate",
                                                  "specific_heating_rate",
                                                  "RT_heating_rate",
                                                  "RT_HI_ionization_rate",
                                                  "RT_HeI_ionization_rate",
                                                  "RT_HeII_ionization_rate",
                                                  "RT_H2_dissociation_rate",
                                                  "H2_self_shielding_length"};

    static_assert(fieldNames.size() == numFields);
    /*struct coolingAttrs
    {
        struct
        {
            int use_grackle;
            int with_radiative_cooling;
            int primordial_chemistry;
            int dust_chemistry;
            int metal_cooling;
            int UVbackground;
            // char* grackle_data_file;
            int    cmb_temperature_floor;
            double Gamma;
            int    h2_on_dust;
            int    use_dust_density_field;
            int    dust_recombination_cooling;
            int    photoelectric_heating;
            double photoelectric_heating_rate;
            int    use_isrf_field;
            double interstellar_radiation_field;
            int    use_volumetric_heating_rate;
            int    use_specific_heating_rate;
            int    three_body_rate;
            int    cie_cooling;
            int    h2_optical_depth_approximation;
            int    ih2co;
            int    ipiht;
            double HydrogenFractionByMass;
            double DeuteriumToHydrogenRatio;
            double SolarMetalFractionByMass;
            double local_dust_to_gas_ratio;
            int    NumberOfTemperatureBins;
            int    CaseBRecombination;
            double TemperatureStart;
            double TemperatureEnd;
            int    NumberOfDustTemperatureBins;
            double DustTemperatureStart;
            double DustTemperatureEnd;
            int    Compton_xray_heating;
            int    LWbackground_sawtooth_suppression;
            double LWbackground_intensity;
            double UVbackground_redshift_on;
            double UVbackground_redshift_off;
            double UVbackground_redshift_fullon;
            double UVbackground_redshift_drop;
            double cloudy_electron_fraction_factor;
            int    use_radiative_transfer;
            int    radiative_transfer_coupled_rate_solver;
            int    radiative_transfer_intermediate_step;
            int    radiative_transfer_hydrogen_only;
            int    self_shielding_method;
            int    H2_self_shielding;
            int    H2_custom_shielding;
            int    h2_charge_exchange_rate;
            int    h2_dust_rate;
            int    h2_h_cooling_rate;
            int    collisional_excitation_rates;
            int    collisional_ionisation_rates;
            int    recombination_cooling_rates;
            int    bremsstrahlung_cooling_rates;
        } attrs;

        auto coolingAttrTuple()
        {
            auto& d = attrs;
            return std::tie(
                d.use_grackle, d.with_radiative_cooling, d.primordial_chemistry, d.dust_chemistry, d.metal_cooling,
                d.UVbackground,
                //!
                //! d.char *grackle_data_file,
                d.cmb_temperature_floor, d.Gamma, d.h2_on_dust, d.use_dust_density_field, d.dust_recombination_cooling,
                d.photoelectric_heating, d.photoelectric_heating_rate, d.use_isrf_field, d.interstellar_radiation_field,
                d.use_volumetric_heating_rate, d.use_specific_heating_rate, d.three_body_rate, d.cie_cooling,
                d.h2_optical_depth_approximation, d.ih2co, d.ipiht, d.HydrogenFractionByMass,
                d.DeuteriumToHydrogenRatio, d.SolarMetalFractionByMass, d.local_dust_to_gas_ratio,
                d.NumberOfTemperatureBins, d.CaseBRecombination, d.TemperatureStart, d.TemperatureEnd,
                d.NumberOfDustTemperatureBins, d.DustTemperatureStart, d.DustTemperatureEnd, d.Compton_xray_heating,
                d.LWbackground_sawtooth_suppression, d.LWbackground_intensity, d.UVbackground_redshift_on,
                d.UVbackground_redshift_off, d.UVbackground_redshift_fullon, d.UVbackground_redshift_drop,
                d.cloudy_electron_fraction_factor, d.use_radiative_transfer, d.radiative_transfer_coupled_rate_solver,
                d.radiative_transfer_intermediate_step, d.radiative_transfer_hydrogen_only, d.self_shielding_method,
                d.H2_self_shielding, d.H2_custom_shielding, d.h2_charge_exchange_rate, d.h2_dust_rate,
                d.h2_h_cooling_rate, d.collisional_excitation_rates, d.collisional_ionisation_rates,
                d.recombination_cooling_rates, d.bremsstrahlung_cooling_rates);
        }

        // Cooling attribute array
        constexpr static std::array attrNames{
            "use_grackle", "with_radiative_cooling", "primordial_chemistry", "dust_chemistry", "metal_cooling",
            "UVbackground",
            //!
            //! d.char *grackle_data_file",
            "cmb_temperature_floor", "Gamma", "h2_on_dust", "use_dust_density_field", "dust_recombination_cooling",
            "photoelectric_heating", "photoelectric_heating_rate", "use_isrf_field", "interstellar_radiation_field",
            "use_volumetric_heating_rate", "use_specific_heating_rate", "three_body_rate", "cie_cooling",
            "h2_optical_depth_approximation", "ih2co", "ipiht", "HydrogenFractionByMass", "DeuteriumToHydrogenRatio",
            "SolarMetalFractionByMass", "local_dust_to_gas_ratio", "NumberOfTemperatureBins", "CaseBRecombination",
            "TemperatureStart", "TemperatureEnd", "NumberOfDustTemperatureBins", "DustTemperatureStart",
            "DustTemperatureEnd", "Compton_xray_heating", "LWbackground_sawtooth_suppression", "LWbackground_intensity",
            "UVbackground_redshift_on", "UVbackground_redshift_off", "UVbackground_redshift_fullon",
            "UVbackground_redshift_drop", "cloudy_electron_fraction_factor", "use_radiative_transfer",
            "radiative_transfer_coupled_rate_solver", "radiative_transfer_intermediate_step",
            "radiative_transfer_hydrogen_only", "self_shielding_method", "H2_self_shielding", "H2_custom_shielding",
            "h2_charge_exchange_rate", "h2_dust_rate", "h2_h_cooling_rate", "collisional_excitation_rates",
            "collisional_ionisation_rates", "recombination_cooling_rates", "bremsstrahlung_cooling_rates"};

        using FieldVariant = std::variant<int*, double*>;

        auto data()
        {
            return std::apply(
                [](auto&... a)
                {
                    auto ret = std::array<FieldVariant, sizeof...(a)>{&a...};
                    std::cout << "ret " << ret.size() << std::endl;
                    return ret;
                },
                coolingAttrTuple());
        }
    } coolingAttrs;*/
    // Units
    T m_code_in_ms  = 1e16;
    T l_code_in_kpc = 46400.;

    //! @brief Unified interface to attribute initialization, reading and writing
    template<class Archive>
    void loadOrStoreAttributes(Archive* ar)
    {
        //! @brief load or store an attribute, skips non-existing attributes on load.
        auto optionalIO = [ar](const std::string& attribute, auto* location, size_t attrSize)
        {
            std::cout << "read " << attribute << std::endl;

            try
            {
                ar->stepAttribute("chem::" + attribute, location, attrSize);
                std::cout << "set " << attribute << std::endl;
            }
            catch (std::out_of_range&)
            {
                std::cout << "Attribute chem::" << attribute << " not set in file, setting to default value "
                          << std::endl;
            }
        };

        optionalIO("m_code_in_ms", &m_code_in_ms, 1);
        optionalIO("l_code_in_kpc", &l_code_in_kpc, 1);
        /*std::apply(
            [&](auto&... attribute)
            {
                std::apply([&](auto&... location) { (optionalIO(std::string(attribute), &location, 1), ...); },
                           coolingAttrs.coolingAttrTuple());
            },
            coolingAttrs.attrNames);*/
    }

private:
    template<size_t... Is>
    auto dataTuple_helper(std::index_sequence<Is...>)
    {
        return std::tie(fields[Is]...);
    }
};

} // namespace cooling
