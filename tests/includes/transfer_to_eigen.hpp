#pragma once

#include <vector>
#include <nlohmann/json.hpp>

#include <Injector/Grids/Defines.h>

using json = nlohmann::json;
using VR = std::vector<GPN::RealType>;

GPN::LogValuesContainer transfer_to_eigen(
    const VR &data, 
    const GPN::RealType factor = 1.0)
{
    using namespace GPN;

    LogValuesContainer out(data.size());
    for (auto i{0ull}; i < data.size(); ++i)
        out(i) = factor * data[i];
    return out;
}