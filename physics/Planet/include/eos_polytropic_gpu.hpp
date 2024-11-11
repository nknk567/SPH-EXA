#pragma once

#include "cstone/sfc/box.hpp"
#include "cstone/traversal/groups.hpp"
#include "cstone/tree/octree.hpp"
#include "cstone/tree/definitions.h"
#include "sph/timestep.h"

namespace planet
{
template<typename T1, typename T2, typename T3, typename Trho, typename Tp, typename Tc>
extern void computePolytropicEOS_HydroStdGPU(size_t firstParticle, size_t lastParticle, T1 Kpoly, T2 exp_poly, T3 gamma,
                                          const Trho* rho, Tp* p, Tc* c);

} // namespace planet
