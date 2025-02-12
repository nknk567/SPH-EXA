//
// Created by Noah Kubli on 11.03.2024.
//

#pragma once

namespace disk
{
template<typename Treal, typename Tmass>
void computeCentralForceGPU(size_t first, size_t last, const Treal* x, const Treal* y, const Treal* z, Treal* ax,
                            Treal* ay, Treal* az, const Tmass* m, StarData& star);
}
