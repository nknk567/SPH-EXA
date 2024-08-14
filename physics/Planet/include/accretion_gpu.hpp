//
// Created by Noah Kubli on 12.03.2024.
//

#pragma once

/*template<typename T1, typename Th, typename Tremove, typename T2>
void computeAccretionConditionGPU(size_t first, size_t last, const T1* x, const T1* y, const T1* z, const Th* h,
                                  Tremove* remove, const T2* spos, T2 star_size, T2 removal_limit_h);*/

//template<typename T1, typename Th, typename Tremove, typename StarData, typename Tm, typename Tv>
template <typename Dataset, typename StarData>
void computeAccretionConditionGPU(size_t first, size_t last, Dataset &d,
                                  //const T1* x, const T1* y, const T1* z, const Th* h,
                                  //Tremove* remove, const Tm* m, const Tv* vx, const Tv* vy, const Tv* vz,
                                  StarData &star);
                                  //const T2* spos, T2 star_size, T2 removal_limit_h, T2& m_accr, T2& vx_accr,
                                  //T2& vy_accr, T2& vz_accr, T2Int &n_removed_local, T2Int &n_accreted_local);

template<typename Tv, typename Tm, typename Tstar>
void sumMassAndMomentumGPU(size_t sum_first, size_t sum_last, const Tv* vx, const Tv* vy, const Tv* vz, const Tm* m,
                           Tv* scratch, Tstar* m_sum, Tstar* p_sum);

template<typename Torder>
void computeNewOrderGPU(size_t first, size_t last, Torder* order, size_t* n_accr, size_t* n_rem2);

template<typename T, typename Torder>
void applyNewOrderGPU(size_t first, size_t last, T* x, T* scratch, Torder* order);
