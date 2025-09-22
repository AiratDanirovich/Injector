#pragma once

#include <vector>
#include <string>

#include <nlohmann/json.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/GridRefiners.hpp>

using json = nlohmann::json;
using VR = std::vector<GPN::RealType>;

const VR make_r_stencils(const json &data, const auto &well_holes)
{
    using namespace GPN;

    const std::string r_grid_type = data["grid"]["r_grid_type"];
    const RealType
        rMin{data["grid"]["r_start"]},
        rMax{data["grid"]["r_end"]};

    if (r_grid_type == "uniform")
    {
        const auto &data2 = data["grid"]["r_uniform_grid"];
        return well_holes.generate_uniform_radial_grid(
            rMin, rMax, data2["rNodes"]);
    }
    else if (r_grid_type == "log")
    {
        const auto &data2 = data["grid"]["r_log_grid"];
        return well_holes.generate_log_radial_grid(
            rMin, rMax, data2["q"].get<RealType>(), data2["r_max_step"].get<RealType>());
    }
    else
        throw std::runtime_error("Incorrect radial grid descriptor.");
}

const GPN::Grids::AbstractRefinerRadial *make_r_refiner(
    const json &data)
{
    using namespace GPN;
    using namespace GPN::Grids;

    const std::string r_grid_type = data["grid"]["r_grid_type"];
    const RealType
        rMin{data["grid"]["r_start"]},
        rMax{data["grid"]["r_end"]};

    if (r_grid_type == "uniform")
    {
        const auto &data2 = data["grid"]["r_uniform_grid"];
        return new RefinerRadial_UniformWellHoles{data2["rNodes"]};
    }
    else if (r_grid_type == "log")
    {
        const auto &data2 = data["grid"]["r_log_grid"];
        return new RefinerRadial_LogWellHoles{
            data2["q"].get<RealType>(),
            data2["r_max_step"].get<RealType>()};
    }
    else
        throw std::runtime_error("Incorrect radial grid descriptor.");
}