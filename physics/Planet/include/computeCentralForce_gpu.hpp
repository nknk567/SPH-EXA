//
// Created by Noah Kubli on 11.03.2024.
//

#pragma once

template<typename Tpos, typename Ta, typename Tm, typename Ts>
void computeCentralForceGPU(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Ta* ax, Ta* ay,
                            Ta* az, const Tm* m, const Ts* spos, Ts sm, Ts* sf, Ts* spot, Tpos g, Ts inner_size);
