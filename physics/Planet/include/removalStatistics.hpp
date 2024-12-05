//
// Created by Noah Kubli on 06.12.2024.
//

#pragma once
#include "cstone/tree/definitions.h"

namespace disk
{
struct RemovalStatistics
{
    double   mass;
    cstone::Vec3<double> momentum;
//    double   momentum[3];
    unsigned count;

    HOST_DEVICE_FUN friend RemovalStatistics operator+(const RemovalStatistics& a, const RemovalStatistics& b)
    {
        RemovalStatistics result;
        result.mass        = a.mass + b.mass;
        result.momentum = a.momentum + b.momentum;
//        result.momentum[0] = a.momentum[0] + b.momentum[0];
//        result.momentum[1] = a.momentum[1] + b.momentum[1];
//        result.momentum[2] = a.momentum[2] + b.momentum[2];
        result.count       = a.count + b.count;
        return result;
    }
};

} // namespace disk