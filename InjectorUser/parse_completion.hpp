#pragma once

#include <array>
#include <limits>

#include <Injector/Grids/Defines.h>

#include <nlohmann/json.hpp>

const std::array<std::array<GPN::RealType, 6>, 6> parse_completion(const nlohmann::json &data)
{
    using namespace GPN;
    using json = nlohmann::json;    

    std::array<std::array<RealType, 6>, 6> out;

    { // flowing fluid
        const auto &data2 = data["fluid"];
        out[0ull] = {data2["density"],
                     data2["specific_heat_capacity"],
                     data2["heat_conductivity"],
                     data["completion"]["tube"]["inner_radius"],
                     0.0,
                     std::numeric_limits<RealType>::max()};
    }

    { // tube
        const auto &data2 = data["completion"]["tube"];
        out[1ull] = {data2["density"],
                     data2["specific_heat_capacity"],
                     data2["heat_conductivity"],
                     data2["thickness"],
                     data2["inner_radius"],
                     data2["depth"]};
    }

    { // annulus
        const auto &data2 = data["completion"]["annulus"];
        out[2ull] = {data2["density"],
                     data2["specific_heat_capacity"],
                     data2["heat_conductivity"],
                     data2["thickness"],
                     // tube inner_radius + tube wall thickness
                     out[1ull][4ull] + out[1ull][3ull],
                     std::numeric_limits<RealType>::max()};
    }

    { // column
        const auto &data2 = data["completion"]["column"];
        out[3ull] = {data2["density"],
                     data2["specific_heat_capacity"],
                     data2["heat_conductivity"],
                     data2["thickness"],
                     // annulus inner_radius + annulus wall thickness
                     out[2ull][4ull] + out[2ull][3ull],
                     std::numeric_limits<RealType>::max()};
    }

    { // cement
        const auto &data2 = data["completion"]["cement"];
        out[4ull] = {data2["density"],
                     data2["specific_heat_capacity"],
                     data2["heat_conductivity"],
                     data2["thickness"],
                     // column inner_radius + column wall thickness
                     out[3ull][4ull] + out[3ull][3ull],
                     std::numeric_limits<RealType>::max()};
    }

    return out;
}