#pragma once

#include <vector>
#include <map>
#include <nlohmann/json.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Completion.hpp>

using json = nlohmann::json;

using VR = std::vector<GPN::RealType>;

const auto get_completion(const json &data, const auto &grid_z)
{
    using namespace std;

    using namespace GPN;
    using namespace GPN::Grids;
    using namespace GPN::Completion;

    map<MaterialType::material_type, VarRing> casing;

    { // column
        const auto &data2 = data["completion"]["column"];
        casing.emplace(
            Completion::MaterialType::Column,
            FactoryVarRing::create_column(
                VarRing{
                    VarRingSimple{
                        VarColumn{
                            Completion::Density{data2["density"].get<VR>()},
                            Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()},
                            Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()}},
                        Completion::VarThickness{data2["radial_thickness"].get<VR>()},
                        Completion::VarDepth{data2["depth_interval"].get<VR>()}},
                    Completion::VarInnerRadius{data2["inner_radius"].get<VR>()}},
                grid_z));
    }
    { // tube
        const auto &data2 = data["completion"]["tube"];
        casing.emplace(
            Completion::MaterialType::Tube,
            FactoryVarRing::create_tube(
                VarRing{
                    VarRingSimple{
                        VarTube{
                            Completion::Density{data2["density"].get<VR>()},
                            Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()},
                            Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()}},
                        Completion::VarThickness{data2["radial_thickness"].get<VR>()},
                        Completion::VarDepth{data2["depth_interval"].get<VR>()}},
                    Completion::VarInnerRadius{data2["inner_radius"].get<VR>()}},
                casing.at(MaterialType::Column),
                grid_z));
    }
    { // flowing fluid
        const auto &data2 = data["fluid"];

        const Flow flow_ring_props{
            GPN::Density{data2["density"].get<RealType>()},
            GPN::SpecificHeatCapacity{data2["specific_heat_capacity"].get<RealType>()},
            GPN::HeatConductivity{data2["heat_conductivity"].get<RealType>()}};

        casing.emplace(
            Completion::MaterialType::Flow,
            FactoryVarRing::create_flow(
                flow_ring_props,
                casing.at(MaterialType::Tube),
                casing.at(MaterialType::Column),
                grid_z));
    }

    { // annulus
        const auto &data2 = data["completion"]["annulus"];

        const VarAnnulus annulus_ring_props{
            Completion::Density{data2["density"].get<VR>()},
            Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()},
            Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()}};

        casing.emplace(
            Completion::MaterialType::Annulus,
            FactoryVarRing::create_annulus(
                annulus_ring_props,
                Completion::VarDepth{data2["depth_interval"].get<VR>()},
                casing.at(MaterialType::Tube),
                casing.at(MaterialType::Column),
                grid_z));
    }

    { // cement
        const auto &data2 = data["completion"]["cement"];
        casing.emplace(Completion::MaterialType::Cement,
                       FactoryVarRing::create_cement(
                           VarRingSimple{
                               VarCement{
                                   Completion::Density{data2["density"].get<VR>()},
                                   Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()},
                                   Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()}},
                               Completion::VarThickness{data2["radial_thickness"].get<VR>()},
                               Completion::VarDepth{data2["depth_interval"].get<VR>()}},
                           casing.at(MaterialType::Column),
                           grid_z));
    }

    return casing;
}