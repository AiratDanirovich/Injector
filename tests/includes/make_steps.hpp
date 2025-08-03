#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>

using VR = std::vector<GPN::RealType>;

VR make_steps(const VR &data)
{
    VR out(data.size() - 1ull);
    for (auto i{1ull}; i < data.size(); ++i)
        out[i - 1ull] = data[i] - data[i - 1ull];

    return out;
}