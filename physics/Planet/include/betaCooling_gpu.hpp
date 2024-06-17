//
// Created by Noah Kubli on 17.04.2024.
//

#pragma once

template<typename Tpos, typename Ts, typename Tdu, typename Trho, typename Tu>
void betaCoolingGPU(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Tdu* du, const Tu* u, Ts star_mass,
                     const Ts* star_pos, Ts beta, Tpos g, Trho *rho, Trho cooling_rho_limit);
