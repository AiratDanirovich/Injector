#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>

using VR = std::vector<GPN::RealType>;

VR transfer_to_vector(
    const Eigen::ArrayX<GPN::RealType> &data)
{
    using namespace GPN;

    VR out(data.size());
    for (auto i{0ll}; i < data.size(); ++i)
        out[i] = data(i);
    return out;
}