//
// Created by Noah Kubli on 17.04.2024.
//

#pragma once

namespace planet
{
template<typename Tpos, typename Tu, typename Ts, typename Tdu, typename Trho, typename Trho2>
void betaCoolingGPU(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, const Tu* u, Tdu* du,
                    Ts star_mass, const Ts* star_pos, Ts beta, Tpos g, const Trho* rho, Ts u_floor,
                    Trho2 cooling_rho_limit);
template<typename Dataset, typename StarData>
void duTimestepGPU(size_t first, size_t last, const Dataset& d, const StarData& star);

} // namespace planet