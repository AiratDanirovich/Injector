/**
 * @file make_history.hpp
 * @brief The file describes the History object creation from a json data
 *
 * 1. By default, an RFP pressure control is assumed. So, it is expected that the
 * nodes "collector->explicit->weights" and "history->rate" exist in the json file.
 * 
 * 2. Bottomhole pressure control is introduced. So, instead of these two nodes one may expect 
 * "history->control_type = bottomhole_pressure" and "history->dynamic->bottomhole_pressure".
 *
 * @author Arthur Salamatin
 * @date 2026-01-21
 */

#pragma once

#include <vector>
#include <string>

#include <nlohmann/json.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/History/History.hpp>

#include "generate_stencils_and_steps.hpp"

using VR = std::vector<GPN::RealType>;

auto make_history(const json &data)
{
    using namespace GPN;
    using namespace GPN::Logs;
#pragma region SET-UNIT-OF-TIME
    const std::string t_unit = data["history"]["t_unit"].get<std::string>();
    RealType factor{1.0};
    if (t_unit == "d")
        factor = 24 * 60 * 60;
    else if (t_unit == "h")
        factor = 60 * 60;
    else if (t_unit == "m")
        factor = 60;
    else if (t_unit == "s")
        factor = 1;
    else
        throw std::runtime_error("Incorrect unit of time.");
#pragma endregion
    const std::string history_type = data["history"]["history_type"];
    if (history_type == "dynamic")
    {
        const auto &data2 = data["history"]["dynamic"];
        VR t_major_steps = data2["t_major_step"].get<VR>();
        for (auto &v : t_major_steps)
            v = v * factor;
        const VR well_rates = data2["well_rate"].get<VR>();
        const VR inlet_temps = data2["inlet_temperature"].get<VR>();

        return HistoryFactory::createFixedRate(t_major_steps, well_rates, inlet_temps);
    }
    else if (history_type == "static")
    {
        const auto &data2 = data["history"]["static"];
        const RealType
            t0{data2["t_start"] * factor},
            t1{data2["t_end"] * factor};

        const RealType t_major_step = std::min(t1 - t0, data2["t_major_step"].get<RealType>() * factor);

        const VR t_stencils{generate_stencils(t0, t1, t_major_step)};

        const RealType well_rate{data2["well_rate"].get<RealType>()}; // m^3/s
        const RealType inlet_temperature{data2["inlet_temperature"].get<RealType>()};

        const std::vector<RealType> t_major_steps{generate_steps(t_stencils)};
        const std::vector<RealType> well_rates(t_major_steps.size(), well_rate);
        const std::vector<RealType> inlet_temperature_set(t_major_steps.size(), inlet_temperature);

        return HistoryFactory::createFixedRate(t_major_steps, well_rates, inlet_temperature_set);
    }
    else
        throw std::runtime_error("Incorrect history type descriptor.");
}