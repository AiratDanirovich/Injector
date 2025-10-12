#pragma once

#include <vector>
#include <exception>

#include <nlohmann/json.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/LogsFactory.hpp>

#include "generate_stencils_and_steps.hpp"

using VR = std::vector<GPN::RealType>;

const GPN::Logs::Geotherma make_geotherma(const json &data, const auto grid2D)
{
    using namespace GPN;

    const auto &data1 = data["collector"]["geotherma"];

    const std::string geotherma_type = data1["type"];
    if (geotherma_type == "const")
    {
        return Logs::GeothermaFactory::create(data1["const"]["initTemperature"], grid2D->first_coord());
    }
    else if (geotherma_type == "interpolate")
    {
        const auto &data2 = data1["interpolate"];
        const VR nodes = data2["z_nodes"].get<VR>();
        const VR vals = data2["t_vals"].get<VR>();
        const RealType z_top = data2["z_top"].get<RealType>();
        return Logs::GeothermaFactory::create(
            nodes,
            vals,
            z_top,
            grid2D->first_coord());
    }
    else
        throw std::runtime_error("Incorrect radial grid descriptors.");
}