#include <fstream>
#include <limits>
#include <iostream>
#include <map>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/Model/Completion.hpp>
#include <Injector/Model/CasingLogs.hpp>
#include <Injector/Model/ExtrudedCasingFactory.hpp>

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/transfer_to_vector.hpp"

using namespace std;

using namespace Catch;
using namespace Catch::Matchers;
using json = nlohmann::json;

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Completion;

using VR = std::vector<RealType>;

TEST_CASE("Well_Test")
{
    ifstream f("heatflow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);

    const auto grid_z{
        Factory::create_axes<CoordinateTypes::Z>(
            RefinerVerticle{
                data["grid"]["z_minor_step"].get<RealType>(),
                transfer_to_eigen(data["collector"]["is_permeable"].get<VR>())},
            Grids::Factory::generate_dual_grid_stencils_from_steps(
                0.0, data["collector"]["thickness"].get<VR>()))};

    map<Completion::MaterialType::material_type, VarRingSimple> casing;

//    casing.emplace(MaterialType::Tube, 1);

    // { // flowing fluid
    //     const auto &data2 = data["fluid"];
    //     casing.push_back(
    //         VarRing{
    //             VarFlow{
    //                 GPN::Density{data2["density"].get<RealType>()},
    //                 GPN::SpecificHeatCapacity{data2["specific_heat_capacity"].get<RealType>()},
    //                 GPN::HeatConductivity{data2["heat_conductivity"].get<RealType>()}},
    //             Completion::VarThickness{data["completion"]["variable"]["tube"]["thickness"].get<VR>()},
    //             Completion::VarDepth{vector{std::numeric_limits<RealType>::max()}}});
    // }
    {
        const auto &data2 = data["completion"]["variable"]["tube"];

        auto density{data2["density"].get<VR>()};
        auto specific_heat_capacity{data2["specific_heat_capacity"].get<VR>()};
        auto heat_conductivity{data2["heat_conductivity"].get<VR>()};
        auto thickness{data2["thickness"].get<VR>()};
    //    auto outer_radius{transfer_to_vector(casing.back().outer_radius)};
        auto depth_interval{data2["depth_interval"].get<VR>()};

        // density.push_back(density.back());
        // specific_heat_capacity.push_back(specific_heat_capacity.back());
        // heat_conductivity.push_back(heat_conductivity.back());
        // thickness.push_back(0.0);
        // outer_radius.push_back(outer_radius.back());
        // depth_interval.push_back(numeric_limits<RealType>::max());

        casing.emplace(Completion::MaterialType::Tube,
            VarRingSimple{
                VarTube{
                    Completion::Density{density},
                    Completion::SpecificHeatCapacity{specific_heat_capacity},
                    Completion::HeatConductivity{heat_conductivity}},
                Completion::VarThickness{thickness},
            //    Completion::VarInnerRadius{outer_radius},
                Completion::VarDepth{depth_interval}
                }
            );
    }
    {
        const auto &data2 = data["completion"]["variable"]["annulus"];
        casing.emplace(Completion::MaterialType::Annulus,
            VarRingSimple{
                VarAnnulus{
                    Completion::Density{data2["density"].get<VR>()},
                    Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()},
                    Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()}},
                Completion::VarThickness{data2["thickness"].get<VR>()},
            //    Completion::VarInnerRadius{casing.back().outer_radius},
                Completion::VarDepth{data2["depth_interval"].get<VR>()}});
    }
    {
        const auto &data2 = data["completion"]["variable"]["column"];
        casing.emplace(Completion::MaterialType::Column,
            VarRingSimple{
                VarColumn{
                    Completion::Density{data2["density"].get<vector<RealType>>()},
                    Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                    Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
                Completion::VarThickness{data2["thickness"].get<vector<RealType>>()},
    //            Completion::VarInnerRadius{casing.back().outer_radius},
                Completion::VarDepth{data2["depth_interval"].get<vector<RealType>>()}});
    }
    {
        const auto &data2 = data["completion"]["variable"]["cement"];
        casing.emplace(Completion::MaterialType::Cement,
            VarRingSimple{
                VarCement{
                    Completion::Density{data2["density"].get<vector<RealType>>()},
                    Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                    Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
                Completion::VarThickness{data2["thickness"].get<vector<RealType>>()},
        //        Completion::VarInnerRadius{casing.back().outer_radius},
                Completion::VarDepth{data2["depth_interval"].get<vector<RealType>>()}});
    }

    // Casing<VarRing> completion{casing};

    // const auto flow_ring =
    //     ExtrudedRingFactory::create_ring(
    //         FlowRing{casing[MaterialType::Flow]},
    //         TubeRing{casing[MaterialType::Tube]},
    //         ColumnRing{casing[MaterialType::Column]},
    //         grid_z);

    // cout << "flow_ring :: thickness:\n"
    //      << flow_ring.thickness.transpose() << endl
    //      << endl;

    // const ExtrudedCasing extr_completion{
    //     ExtrudedCasingFactory::create(completion, grid_z)};
}
