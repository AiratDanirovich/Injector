#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>

using VR = std::vector<GPN::RealType>;


GPN::LogValuesContainer set_is_permeable_stencils(
    const GPN::LogValuesContainer &perforated_data, 
    const std::vector<std::ptrdiff_t> &ids)
{
    using namespace GPN;

    LogValuesContainer out(perforated_data);
    for (const auto i : ids)
    {    
        out(i) = 1.0;
    }

    return out;
}