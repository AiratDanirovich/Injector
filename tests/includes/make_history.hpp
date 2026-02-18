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

auto read_units_factor(const json &data)
{
    using namespace GPN;
    const std::string t_unit{data["history"]["t_unit"].get<std::string>()};
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

    return factor;
}

auto read_z_ref(const json &data)
{
    using namespace GPN;
    try
    {
        return data["collector"]["hydrostat"]["z_ref"].get<RealType>();
    }
    catch(...)
    {
        return 100.0;
    }
}

auto make_history(const json &data)
{
    using namespace GPN;
    using namespace GPN::Logs;
    
    const RealType Atm2Pa{1e5};
    const RealType factor{read_units_factor(data)};
#pragma region CHOOSE-HISTORY-TYPE
    const std::string history_type{data["history"]["history_type"].get<std::string>()};
    const std::string control_type{data["history"].value<std::string>("control_type", "RFP")};
    if (history_type == "dynamic")
    {
        const auto &data2 {data["history"]["dynamic"]};
        VR t_major_steps{data2["t_major_step"].get<VR>()};
        for (auto &v : t_major_steps)
            v = v * factor; // change units of time-steps to seconds
        const VR inlet_temps {data2["inlet_temperature"].get<VR>()};
        if ((control_type == "RFP") || (control_type == "bottomhole_rate"))
        {
            const VR well_rates { data2["well_rate"].get<VR>()};
            return HistoryFactory::createFixedRate(t_major_steps, well_rates, inlet_temps);
        }
        else if (control_type == "bottomhole_pressure")
        {
            VR bothole_pressure = data2["bottomhole_pressure"].get<VR>();
            for (auto &v : bothole_pressure)
                v = v * Atm2Pa; // change units of pressure to Pa
            return HistoryFactory::createFixedPressure(t_major_steps, bothole_pressure, inlet_temps, read_z_ref(data));
        }
        else
            throw std::runtime_error("Incorrect history control type.");
    }
    else if (history_type == "static")
    {
        const auto &data2 = data["history"]["static"];
        const RealType
            t0{data2["t_start"] * factor},
            t1{data2["t_end"] * factor};

        const RealType t_major_step = std::min(t1 - t0, data2["t_major_step"].get<RealType>() * factor);

        const VR t_stencils{generate_stencils(t0, t1, t_major_step)};
        const std::vector<RealType> t_major_steps{generate_steps(t_stencils)};

        const RealType inlet_temperature{data2["inlet_temperature"].get<RealType>()};
        const std::vector<RealType> inlet_temps(t_major_steps.size(), inlet_temperature);

        if ((control_type == "RFP") || (control_type == "bottomhole_rate"))
        {
            const RealType well_rate{data2["well_rate"].get<RealType>()}; // m^3/s
            const std::vector<RealType> well_rates(t_major_steps.size(), well_rate);
            return HistoryFactory::createFixedRate(t_major_steps, well_rates, inlet_temps);
        }
        else if (control_type == "bottomhole_pressure")
        {
            const RealType bothole_pressure{Atm2Pa*data2["bothole_pressure"].get<RealType>()}; // 
            const std::vector<RealType> bothole_pressures(t_major_steps.size(), bothole_pressure);
            return HistoryFactory::createFixedPressure(t_major_steps, bothole_pressures, inlet_temps, read_z_ref(data));
        }
        else
            throw std::runtime_error("Incorrect history control type.");
    }
    else
        throw std::runtime_error("Incorrect history type descriptor.");
#pragma endregion
}

auto read_minor_step(const json &data)
{
    using namespace GPN;
    return data["history"]["t_minor_step"].get<RealType>()*read_units_factor(data);
}

auto read_start_time(const json &data)
{
    using namespace GPN;
    return data["history"]["start_time"].get<RealType>()*read_units_factor(data);
}