//
// Created by Noah Kubli on 17.04.2024.
//

#pragma once

namespace disk
{

template<typename Treal, typename Thydro, typename StarData>
void betaCoolingGPU(size_t first, size_t last, const Treal* x, const Treal* y, const Treal* z, const Treal* u,
                    const Thydro* rho, Treal* du, StarData& star);

template<typename Treal, typename StarData>
double duTimestepGPU(size_t first, size_t last, const Treal* u, const Treal* du);

} // namespace disk
