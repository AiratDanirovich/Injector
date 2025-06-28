#pragma once

#include <vector>
#include <nlohmann/json.hpp>

#include <Injector/Model/Completion.hpp>

using json = nlohmann::json;

const auto get_completion(const json &data)
{
    using namespace GPN;
    using namespace GPN::Completion;

    std::vector<Ring> out;
    out.reserve(6);

    { // flowing fluid
        const auto &data2 = data["fluid"];
        out.push_back(
            Ring{
                Flow{
                    Density{data2["density"]},
                    SpecificHeatCapacity{data2["specific_heat_capacity"]},
                    GPN::HeatConductivity{data2["heat_conductivity"]}},
                Thickness{data["completion"]["tube"]["inner_radius"]},
                InnerRadius{0.0},
                Depth{std::numeric_limits<RealType>::max()}});
    }

    { // tube
        const auto &data2 = data["completion"]["tube"];

        out.push_back(
            Ring{
                Tube{
                    Density{data2["density"]},
                    SpecificHeatCapacity{data2["specific_heat_capacity"]},
                    GPN::HeatConductivity{data2["heat_conductivity"]}},
                Thickness{data2["thickness"]},
                InnerRadius{out.back().outer_radius},
                Depth{data2["depth"]}});
    }

    { // annulus
        const auto &data2 = data["completion"]["annulus"];

        out.push_back(
            Ring{
                Annulus{
                    Density{data2["density"]},
                    SpecificHeatCapacity{data2["specific_heat_capacity"]},
                    GPN::HeatConductivity{data2["heat_conductivity"]}},
                Thickness{data2["thickness"]},
                InnerRadius{out.back().outer_radius},
                Depth{data2["depth"]}});
    }

    { // column
        const auto &data2 = data["completion"]["column"];

        out.push_back(
            Ring{
                Annulus{
                    Density{data2["density"]},
                    SpecificHeatCapacity{data2["specific_heat_capacity"]},
                    GPN::HeatConductivity{data2["heat_conductivity"]}},
                Thickness{data2["thickness"]},
                InnerRadius{out.back().outer_radius},
                Depth{data2["depth"]}});
    }

    { // cement_1
        const auto &data2 = data["completion"]["cement"]["inner"];

        out.push_back(
            Ring{
                Cement{
                    Density{data2["density"]},
                    SpecificHeatCapacity{data2["specific_heat_capacity"]},
                    GPN::HeatConductivity{data2["heat_conductivity"]}},
                Thickness{data2["thickness"]},
                InnerRadius{out.back().outer_radius},
                Depth{data2["depth"]}});
    }

    { // cement_2
        const auto &data2 = data["completion"]["cement"]["outer"];

        out.push_back(
            Ring{
                Cement{
                    Density{data2["density"]},
                    SpecificHeatCapacity{data2["specific_heat_capacity"]},
                    GPN::HeatConductivity{data2["heat_conductivity"]}},
                Thickness{data2["thickness"]},
                InnerRadius{out.back().outer_radius},
                Depth{data2["depth"]}});
    }

    return out;
}