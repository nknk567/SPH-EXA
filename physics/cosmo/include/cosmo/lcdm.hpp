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
 * @brief Lambda Cold Dark Matter (CDM) Cosmology.
 *
 * @author Jonathan Coles <jonathan.coles@cscs.ch>
 * @author Noah Kubli <noah.kubli@uzh.ch>
 */

#pragma once

#include <cmath>
#include "utils.hpp"

namespace cosmo
{

template<typename T>
struct LambdaCDMParameters
{
    //! @brief Hubble Constant today [km/Mpc/s]
    T H0{0.0};
    //! @brief Matter Density
    T OmegaMatter{0.0};
    //! @brief Radiation Density
    T OmegaRadiation{0.0};
    //! @brief Cosmology Constant (vacuum density)
    T OmegaLambda{0.0};

    //        T H0, OmegaMatter, OmegaRadiation, OmegaLambda;
    inline constexpr static std::array names{"H0", "OmegaMatter", "OmegaRadiation", "OmegaLambda"};
    // Assert T is in FieldVariant
    auto dataTuple()
    {
        auto ret = std::tie(H0, OmegaMatter, OmegaRadiation, OmegaLambda);

        static_assert(std::tuple_size_v<decltype(ret)> == names.size());
        return ret;
    }
};

static constexpr LambdaCDMParameters<double> Planck2018 = {
    .H0             = 67.66,
    .OmegaMatter    = 0.3111,
    .OmegaRadiation = 0.0,
    .OmegaLambda    = 0.6889,
};

template<typename T, int rel_prec = -7, int abs_prec = -7>
class LambdaCDM : public Cosmology<T>
{
public:
    T relativeError = std::pow(10, rel_prec);
    T absoluteError = std::pow(10, abs_prec);

private:
    using FieldVariant = Cosmology<T>::FieldVariant;
    using Parameters   = LambdaCDMParameters<T>;

    Parameters parameters = Planck2018;

    std::vector<FieldVariant> getParameters() override
    {
        return std::apply([](auto&... a) { return std::vector<FieldVariant>{&a...}; }, parameters.dataTuple());
    }

    std::vector<const char*> getParameterNames() override
    {
        auto a =
            std::apply([](auto&... a) { return std::array<const char*, sizeof...(a)>{(&a[0])...}; }, Parameters::names);
        return std::vector(a.begin(), a.end());
    }

    static constexpr int cosmology_type = 0;
    int                  getCosmologyType() override { return cosmology_type; }

public:
    // maybe
    // void sanityCheckParameters()

    //
    //! @brief Values from 2018 Planck final data release, TT,TE,EE+lowE+lensing+BAO
    //

    //    LambdaCDM() = default;

    //    LambdaCDM(struct Parameters p)
    //        : LambdaCDM(p.H0, p.OmegaMatter, p.OmegaRadiation, p.OmegaLambda)
    //    {
    //    }
    //
    //    LambdaCDM(T H0, T OmegaMatter, T OmegaRadiation, T OmegaLambda)
    //        : H0(H0)
    //        , OmegaMatter(OmegaMatter)
    //        , OmegaRadiation(OmegaRadiation)
    //        , OmegaLambda(OmegaLambda)
    //    {
    //        if (H0 <= 0) throw std::domain_error("Hubble0 parameter must be strictly positive.");
    //
    //        if (OmegaMatter < 0) throw std::domain_error("OmegaMatter parameter must be positive.");
    //
    //        if (OmegaRadiation < 0) throw std::domain_error("OmegaRadiation parameter must be positive.");
    //
    //        if (OmegaLambda < 0) throw std::domain_error("OmegaLambda parameter must be positive.");
    //
    //        if (OmegaRadiation == 0 && OmegaLambda == 0)
    //            throw std::domain_error(
    //                "In LambdaCDM at least one of OmegaRadiation or OmegaLambda must be positive. Otherwise use
    //                CDM.");
    //
    //        if (OmegaMatter == 0)
    //        {
    //            if (OmegaRadiation == 0)
    //                throw std::domain_error("OmegaRadiation parameter must be strictly positive when OmegaMatter is
    //                zero.");
    //        }
    //    }

    //  template<class Archive>
    //  void loadOrStoreAttributes(Archive* ar)
    //  {
    //      //! @brief load or store an attribute, skips non-existing attributes on load.
    //      auto optionalIO = [ar](const std::string& attribute, auto* location, size_t attrSize)
    //      {
    //          try
    //          {
    //              ar->stepAttribute(attribute, location, attrSize);
    //          }
    //          catch (std::out_of_range&)
    //          {
    //              std::cout << "Attribute " << attribute << " not set in file, setting to default value " << *location
    //                        << std::endl;
    //          }
    //      };

    //      optionalIO("Hubble0", &H0, 1);
    //      optionalIO("OmegaMatter", &OmegaMatter, 1);
    //      optionalIO("OmegaRadiation", &OmegaRadiation, 1);
    //      optionalIO("OmegaLambda", &OmegaLambda, 1);
    //  }

    T hubble_H(const T a) const
    {
        T Omega0         = parameters.OmegaMatter + parameters.OmegaRadiation + parameters.OmegaLambda;
        T OmegaCurvature = T(1) - Omega0;

        T a2 = a * a;
        T a3 = a * a2;
        T a4 = a * a3;
        return (parameters.H0 / a2) * std::sqrt(parameters.OmegaRadiation + parameters.OmegaMatter * a +
                                                OmegaCurvature * a2 + parameters.OmegaLambda * a4);
    }

    T time(const T a) const
    {
        auto f = [this](const T x) { return 2. / (3. * x * hubble_H(pow(x, 2. / 3.))); };
        return integrate_romberg<T>(f, 0.0, std::pow(a, 1.5), 1e-2 * relativeError);
    }

    T scale_factor_a(T t) const
    {
        auto f  = [this, t](const T a) { return time(a) - t; };
        auto df = [this](const T a) { return 1.0 / (a * hubble_H(a)); };

        return find_zero_newton(f, df, t * parameters.H0, 0.0, 1.0e38, relativeError, absoluteError);
    }

    T driftTimeCorrection(T t, T dt) override
    {
        auto f = [this](const T x) { return -x / hubble_H(T(1.0) / x); };
        return integrate_romberg<T>(f, T(1.0) / scale_factor_a(t), T(1.0) / scale_factor_a(t + dt), relativeError);
    }

    T kickTimeCorrection(T t, T dt) override
    {
        auto integrand_transformed = [this](const T x) { return -T(1.0) / hubble_H(T(1.0) / x); };

        auto    trans_f = [](T x) { return T(1.0) / x; };
        const T lower   = trans_f(scale_factor_a(t));
        const T upper   = trans_f(scale_factor_a(t + dt));

        return integrate_romberg<T>(integrand_transformed, lower, upper, relativeError);
    }

    //! @brief Assuming gamma > 1
    T kickTimeCorrectionSPH(T t, T dt, T gamma) const override
    {
        if (gamma <= T(1.0)) throw std::runtime_error("gamma has to be > 1.");

        // We integrate 1. / (a^(3*gamma - 2) * H) da
        // The integrand is transformed in order to remove the inverse a dependence (for numerical stability)
        const T beta                  = T(3.) * gamma - T(2.);
        const T zeta                  = T(1.) - beta;
        auto    integrand_transformed = [this, zeta](T x) { return -T(1.0) / hubble_H(std::pow(x, zeta) * zeta); };

        auto    trans_f = [zeta](T x) { return zeta * std::pow(x, zeta); };
        const T lower   = trans_f(scale_factor_a(t));
        const T upper   = trans_f(scale_factor_a(t + dt));

        return integrate_romberg<T>(integrand_transformed, lower, upper, relativeError);
    }
};

} // namespace cosmo
