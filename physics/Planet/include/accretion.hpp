//
// Created by Noah Kubli on 15.03.2024.
//

#pragma once

#include <cassert>
#include <mpi.h>
#include "cstone/primitives/mpi_wrappers.hpp"
#include "cstone/fields/field_get.hpp"
#include "cstone/tree/accel_switch.hpp"
#include "cstone/cuda/cuda_stubs.h"
#include "cstone/util/type_list.hpp"

#include "sph/particles_data.hpp"

#include "accretion_impl.hpp"
#include "accretion_gpu.hpp"
#include "fieldListExclude.hpp"

namespace planet
{

//! @brief Flag particles for removal. Overwrites keys.
template<typename Dataset, typename StarData>
void computeAccretionCondition(size_t first, size_t last, Dataset& d, StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        computeAccretionConditionGPU(first, last, d,
                                     //rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z),
                                     //rawPtr(d.devData.h), rawPtr(d.devData.keys), rawPtr(d.devData.m),
                                     //rawPtr(d.devData.vx), rawPtr(d.devData.vy), rawPtr(d.devData.vz),
                                     star);
                                     /*star.position.data(), star.inner_size,
                                     star.removal_limit_h, star.m_accreted_local, star.p_accreted_local[0],
                                     star.p_accreted_local[1],
                                     star.p_accreted_local[2], star.n_removed_local,
                                     star.n_accreted_local);*/

    }
    else
    {
        computeAccretionConditionImpl(first, last, d.x.data(), d.y.data(), d.z.data(), d.h.data(), d.keys.data(),
                                      star.position.data(), star.inner_size, star.removal_limit_h);
    }
}

//! @brief Compute new particle ordering with the particles to remove at the end. Overwrites keys.
template<typename Dataset, typename StarData>
void computeNewOrder(size_t first, size_t last, Dataset& d, StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        computeNewOrderGPU(first, last, rawPtr(get<"keys">(d)), &star.n_accreted_local, &star.n_removed_local);
    }
    else { computeNewOrderImpl(first, last, d.keys.data(), &star.n_accreted_local, &star.n_removed_local); }
}

//! @brief Apply the new particle ordering to all the conserved fields. Dependent fields are used as scratch buffer. The
//! particles are ordered as follows: active particles, particles to be accreted, other particles to be removed
template<typename ConservedFields, typename DependentFields, typename Dataset>
void applyNewOrder(size_t first, size_t last, Dataset& d)
{
    using SortFields = decltype(util::FieldList<"x", "y", "z", "h", "m">{} + ConservedFields{});
    auto sortVectors = get<SortFields>(d);

    using ScratchFields       = FieldListExclude_t<"keys", DependentFields>;
    auto scratch_fields_tuple = [&d]()
    {
        if constexpr (std::tuple_size_v<decltype(util::make_array(ScratchFields{}))> == 1)
        {
            return std::tie(get<ScratchFields>(d));
        }
        else { return get<ScratchFields>(d); }
    };
    auto scratchVectors = scratch_fields_tuple();

    util::for_each_tuple(
        [&](auto& vector)
        {
            auto& scratch = util::pickType<std::decay_t<decltype(vector)>&>(scratchVectors);
            if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
            {
                applyNewOrderGPU(first, last, rawPtr(vector), rawPtr(scratch), rawPtr(get<"keys">(d)));
            }
            else { applyNewOrderImpl(first, last, vector.data(), scratch.data(), d.keys.data()); }
        },
        sortVectors);
}

//! @brief Compute total angular momentum and mass of the accreted particles. Dependent Fields are used as scratch
//! buffers.
template<typename DependentFields, typename Dataset, typename StarData>
void sumAccretedMassAndMomentum(size_t first, size_t last, Dataset& d, StarData& star)
{
    if (star.n_accreted_local + star.n_removed_local > (last - first))
        throw std::runtime_error("Accreting more particles than on rank");
    auto  scratchBuffers = get<DependentFields>(d);
    auto& scratch        = util::pickType<std::decay_t<decltype(get<"vx">(d))>&>(scratchBuffers);

    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        sumMassAndMomentumGPU(last - star.n_accreted_local - star.n_removed_local, last - star.n_removed_local,
                              rawPtr(get<"vx">(d)), rawPtr(get<"vy">(d)), rawPtr(get<"vz">(d)), rawPtr(get<"m">(d)),
                              rawPtr(scratch), &star.m_accreted_local, star.p_accreted_local.data());
    }
    else
    {
        sumMassAndMomentumImpl(last - star.n_accreted_local - star.n_removed_local, last - star.n_removed_local,
                               get<"vx">(d).data(), get<"vy">(d).data(), get<"vz">(d).data(), get<"m">(d).data(),
                               scratch.data(), &star.m_accreted_local, star.p_accreted_local.data());
    }
}

//! @brief Exchange accreted mass and momentum between ranks and add to star.
template<typename StarData>
void exchangeAndAccreteOnStar(StarData& star, double minDt_m1, int rank)
{
    double                m_accreted_global{};
    std::array<double, 3> p_accreted_global{};

    MPI_Reduce(&star.m_accreted_local, &m_accreted_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(star.p_accreted_local.data(), p_accreted_global.data(), 3, MpiType<double>{}, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (rank == 0)
    {
        double m_star_new = m_accreted_global + star.m;

        std::array<double, 3> p_star;
        for (size_t i = 0; i < 3; i++)
        {
            p_star[i] = (star.position_m1[i] / minDt_m1) * star.m;
            p_star[i] += p_accreted_global[i];
            star.position_m1[i] = p_star[i] / m_star_new * minDt_m1;
        }

        star.m = m_star_new;

        printf("accreted mass: %g\tstar mass: %g\n", m_accreted_global, star.m);
        printf("accreted momentum x: %g\tstar momentum x: %g\n", p_accreted_global[0], p_star[0]);
        printf("accreted momentum y: %g\tstar momentum y: %g\n", p_accreted_global[1], p_star[1]);
        printf("accreted momentum z: %g\tstar momentum z: %g\n", p_accreted_global[2], p_star[2]);
    }

    if (rank == 0) { printf("accreted mass local: %g\n", star.m_accreted_local); }

    MPI_Bcast(star.position_m1.data(), 3, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&star.m, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
}

} // namespace planet
