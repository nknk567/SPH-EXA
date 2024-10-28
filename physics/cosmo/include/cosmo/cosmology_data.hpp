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
 * @brief Cosmological parameters and associated functions.
 *
 * @author Jonathan Coles <jonathan.coles@cscs.ch>
 */

#pragma once

namespace cosmo
{

template<typename T>
class Cosmology // : public cstone::FieldStates<CosmologyData<T>>
{
public:
    virtual T   driftTimeCorrection(T t, T dt)                  = 0;
    virtual T   kickTimeCorrection(T t, T dt)                   = 0;
    virtual T   kickTimeCorrectionSPH(T t, T dt, T gamma) const = 0;
    virtual int getCosmologyType()                              = 0;
    virtual ~Cosmology()                                        = default;

    template<typename Archive>
    void loadOrStoreAttributes(Archive* ar)
    {
        //! @brief load or store an attribute, skips non-existing attributes on load.
        auto optionalIO = [ar](const std::string& attribute, auto* location, size_t attrSize)
        {
            try
            {
                ar->stepAttribute("cosmology::" + attribute, location, attrSize);
            }
            catch (std::out_of_range&)
            {
                if (ar->rank() == 0)
                {
                    std::cout << "Attribute cosmology::" << attribute
                              << " not set in file or initializer, setting to default value " << *location << std::endl;
                }
            }
        };

        // This is a readonly attribute here, it is set in the factory method.
        int cosmology_type = getCosmologyType();
        optionalIO("cosmology_type", &cosmology_type, 1);

        auto parameterNames = getParameterNames();
        auto parameters     = getParameters();
        for (size_t i = 0; i < parameterNames.size(); i++)
        {
            std::visit([&](auto* location) { optionalIO(std::string(parameterNames[i]), location, 1); }, parameters[i]);
        }
    }

protected:
    using FieldVariant = std::variant<double*, float*>;

private:
    virtual std::vector<const char*>  getParameterNames();
    virtual std::vector<FieldVariant> getParameters();
};

} // namespace cosmo
