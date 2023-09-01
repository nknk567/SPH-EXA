//
// Created by Noah Kubli on 09.07.23.
//

#pragma once

namespace cooling
{

//! @brief the maximum time-step based on local particles that Grackle can tolerate
template<class Dataset, typename Cooler, typename Chem>
auto coolingTimestep(size_t first, size_t last, Dataset& d, Cooler& cooler, Chem& chem)
{
    using T             = typename Dataset::RealType;
    using CoolingFields = typename Cooler::CoolingFields;

    T minTc(INFINITY);
#pragma omp parallel for reduction(min : minTc)
    for (size_t i = first; i < last; i++)
    {
        const T cooling_time = cooler.cooling_time(d.rho[i], d.u[i], cstone::getPointers(get<CoolingFields>(chem), i));
        d.ct[i]              = cooling_time;
        minTc                = std::min(std::abs(cooler.ct_crit * cooling_time), minTc);
    }
    return minTc;
}

template<typename HydroData, typename ChemData, typename Cooler>
void eos_cooling(size_t startIndex, size_t endIndex, HydroData& d, ChemData& chem, Cooler& cooler)
{
    using CoolingFields = typename Cooler::CoolingFields;
    using T             = typename HydroData::RealType;
    const auto* rho     = d.rho.data();

    auto* p = d.p.data();
    auto* c = d.c.data();

#pragma omp parallel for schedule(static)
    for (size_t i = startIndex; i < endIndex; ++i)
    {
        T pressure    = cooler.pressure(rho[i], d.u[i], cstone::getPointers(get<CoolingFields>(chem), i));
        T gamma       = cooler.adiabatic_index(rho[i], d.u[i], cstone::getPointers(get<CoolingFields>(chem), i));
        T sound_speed = std::sqrt(gamma * pressure / rho[i]);
        p[i]          = pressure;
        c[i]          = sound_speed;
    }
}
template<typename HydroData, typename ChemData, typename Cooler>
void eos_cooling_ve(size_t startIndex, size_t endIndex, HydroData& d, ChemData& chem, Cooler& cooler)
{
    using CoolingFields = typename Cooler::CoolingFields;
    using T             = typename HydroData::RealType;

    const auto* u     = d.u.data();
    const auto* m     = d.m.data();
    const auto* kx    = d.kx.data();
    const auto* xm    = d.xm.data();
    const auto* gradh = d.gradh.data();

    auto* prho = d.prho.data();
    auto* c    = d.c.data();

    bool storeRho = (d.rho.size() == d.m.size());
    bool storeP   = (d.p.size() == d.m.size());

    auto* p   = d.p.data();
    auto* rho = d.rho.data();

#pragma omp parallel for schedule(static)
    for (size_t i = startIndex; i < endIndex; ++i)
    {
        auto rho = kx[i] * m[i] / xm[i];
        if (storeRho) { d.rho[i] = rho; }

        T pressure = cooler.pressure(rho, d.u[i], cstone::getPointers(get<CoolingFields>(chem), i));
        if (storeP) { d.p[i] = pressure; }

        T gamma       = cooler.adiabatic_index(rho, d.u[i], cstone::getPointers(get<CoolingFields>(chem), i));
        T sound_speed = std::sqrt(gamma * pressure / rho);
        c[i]          = sound_speed;
        prho[i]       = pressure / (kx[i] * m[i] * m[i] * gradh[i]);
    }
}
} // namespace cooling