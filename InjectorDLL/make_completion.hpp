#pragma once

#include <vector>
#include <array>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Completion.hpp>

const auto make_completion(const std::array<std::array<GPN::RealType, 6>, 6> &data)
{
    using namespace GPN;
    using namespace GPN::Completion;

    std::vector<Ring> out;
    out.reserve(6);

    { // flowing fluid
        const auto &data2 = data[0ull];
        out.push_back(
            Ring{
                Flow{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{std::numeric_limits<RealType>::max()}});
    }

    { // tube
        const auto &data2 = data[1ull];
        out.push_back(
            Ring{
                Tube{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{data2[5ull]}});
    }

    { // annulus
        const auto &data2 = data[2ull];
        out.push_back(
            Ring{
                Annulus{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{std::numeric_limits<RealType>::max()}});
    }

    { // column
        const auto &data2 = data[3ull];
        out.push_back(
            Ring{
                Column{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{std::numeric_limits<RealType>::max()}});
    }

    { // cement_1
        const auto &data2 = data[4ull];
        out.push_back(
            Ring{
                Cement{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{std::numeric_limits<RealType>::max()}});
    }

    return out;
}