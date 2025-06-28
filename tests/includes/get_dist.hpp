#pragma once

#include <algorithm>

#include <Injector/Grids/Defines.h>

std::ptrdiff_t get_dist(const auto &data, GPN::RealType val)
{
    const auto it =
        std::upper_bound(
            data.cbegin(),
            data.cend(),
            val);
    return std::distance(data.cbegin(), it) - 1ll;
}