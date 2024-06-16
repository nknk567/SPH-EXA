//
// Created by Noah Kubli on 12.03.2024.
//

#pragma once

template<typename T1, typename Th, typename Tremove, typename T2>
void computeAccretionConditionGPU(size_t first, size_t last, const T1* x, const T1* y, const T1* z, const Th* h,
                                  Tremove* remove, const T2* spos, T2 star_size, T2 removal_limit_h);

template<typename Tv, typename Tm, typename Tstar>
void sumMassAndMomentumGPU(size_t sum_first, size_t sum_last, const Tv* vx, const Tv* vy, const Tv* vz, const Tm* m,
                           Tv* scratch, Tstar* m_sum, Tstar* p_sum);

template<typename Torder>
void computeNewOrderGPU(size_t first, size_t last, Torder* order, size_t* n_accr, size_t* n_rem2);

template<typename T, typename Torder>
void applyNewOrderGPU(size_t first, size_t last, T* x, T* scratch, Torder* order);
