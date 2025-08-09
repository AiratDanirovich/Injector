#include <fstream>
#include <limits>
#include <iostream>
#include <map>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/Model/Completion.hpp>
#include <Injector/Model/ExtrudedCasingFactory.hpp>

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/transfer_to_vector.hpp"
#include "includes/accumulate_steps.hpp"

using namespace std;

using namespace Catch;
using namespace Catch::Matchers;
using json = nlohmann::json;

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Completion;

using VR = vector<RealType>;

const RealType tol{1e-12};

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

    map<Completion::MaterialType::material_type, VarRing> casing;

    {
        const auto &data2 = data["completion"]["variable"]["column"];
        casing.emplace(
            Completion::MaterialType::Column,
            FactoryVarRing::create_column_ring(
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

        const auto density = Completion::Density{data2["density"].get<VR>()};
        const auto capacity = Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()};
        const auto conductivity = Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()};
        const auto thickness = Completion::VarThickness{data2["radial_thickness"].get<VR>()};
        const auto depth_interval = Completion::VarDepth{data2["depth_interval"].get<VR>()};
        const auto depth_stencils = accumulate_steps(depth_interval.value.data);
        const auto inner_radius = Completion::VarInnerRadius{data2["inner_radius"].get<VR>()};
        const auto ring{casing.at(MaterialType::Column)};
        auto depth_id{1ll};
        for (auto id{0ll}; id < grid_z.mesh_size(); ++id)
        {
            if (grid_z.mesh_nodes(id) > depth_stencils[depth_id])
                ++depth_id;
            REQUIRE(grid_z.mesh_nodes(id) < depth_stencils[depth_id]);

            CHECK(ring.density()(id) == density.value(depth_id - 1ll));
            CHECK(ring.specific_heat_capacity()(id) == capacity.value(depth_id - 1ll));
            CHECK(ring.heat_conductivity()(id) == conductivity.value(depth_id - 1ll));

            CHECK(ring.inner_radius(id) == inner_radius.value(depth_id - 1ll));
            CHECK_THAT(ring.thickness(id), WithinRel(thickness.value(depth_id - 1ll), tol));
            CHECK(ring.inner_radius(id) + ring.thickness(id) == ring.outer_radius(id));
        }
        REQUIRE(depth_id == depth_stencils.size() - 1ull);
        REQUIRE(depth_id == density.value.size());
    }

    {
        const auto &data2 = data["completion"]["variable"]["tube"];
        casing.emplace(
            Completion::MaterialType::Tube,
            FactoryVarRing::create_tube_ring(
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

        const auto density = Completion::Density{data2["density"].get<VR>()};
        const auto capacity = Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()};
        const auto conductivity = Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()};
        const auto thickness = Completion::VarThickness{data2["radial_thickness"].get<VR>()};
        const auto depth_interval = Completion::VarDepth{data2["depth_interval"].get<VR>()};
        const auto depth_stencils = accumulate_steps(depth_interval.value.data);
        const auto inner_radius = Completion::VarInnerRadius{data2["inner_radius"].get<VR>()};

        const auto ring{casing.at(MaterialType::Tube)};
        const auto &column{casing.at(MaterialType::Column)};
        auto id{0ll}, depth_id{1ll};
        for (;
             (id < grid_z.mesh_size()) &&
             (grid_z.mesh_nodes(id) < depth_stencils.back());
             ++id)
        {
            if (grid_z.mesh_nodes(id) > depth_stencils[depth_id])
                ++depth_id;
            REQUIRE(grid_z.mesh_nodes(id) < depth_stencils[depth_id]);

            CHECK(ring.density()(id) == density.value(depth_id - 1ll));
            CHECK(ring.specific_heat_capacity()(id) == capacity.value(depth_id - 1ll));
            CHECK(ring.heat_conductivity()(id) == conductivity.value(depth_id - 1ll));

            CHECK(ring.inner_radius(id) == inner_radius.value(depth_id - 1ll));
            CHECK(ring.thickness(id) == thickness.value(depth_id - 1ll));
            CHECK(ring.inner_radius(id) + ring.thickness(id) == ring.outer_radius(id));
        }
        REQUIRE(depth_id == depth_stencils.size() - 1ull);
        REQUIRE(depth_id == density.value.size());
        for (;
             (id < grid_z.mesh_size()) &&
             (grid_z.mesh_nodes(id) < depth_stencils.back());
             ++id)
        {
            CHECK(ring.density()(id) == density.value(depth_id - 1ll));
            CHECK(ring.specific_heat_capacity()(id) == capacity.value(depth_id - 1ll));
            CHECK(ring.heat_conductivity()(id) == conductivity.value(depth_id - 1ll));

            CHECK(ring.inner_radius(id) == column.inner_radius(id));
            CHECK(ring.outer_radius(id) == column.inner_radius(id));
            CHECK(ring.thickness(id) == 0.0);
        }
    }

    { // flowing fluid
        const auto &data2 = data["fluid"];

        const Flow flow_ring_props{
            GPN::Density{data2["density"].get<RealType>()},
            GPN::SpecificHeatCapacity{data2["specific_heat_capacity"].get<RealType>()},
            GPN::HeatConductivity{data2["heat_conductivity"].get<RealType>()}};

        casing.emplace(
            Completion::MaterialType::Flow,
            FactoryVarRing::create_flow_ring(
                flow_ring_props,
                casing.at(MaterialType::Tube),
                casing.at(MaterialType::Column),
                grid_z));

        const auto density = GPN::Density{data2["density"].get<RealType>()};
        const auto capacity = GPN::SpecificHeatCapacity{data2["specific_heat_capacity"].get<RealType>()};
        const auto conductivity = GPN::HeatConductivity{data2["heat_conductivity"].get<RealType>()};

        const auto ring{casing.at(MaterialType::Flow)};
        const auto &tube{casing.at(MaterialType::Tube)};
        const auto &column{casing.at(MaterialType::Column)};
        auto id{0ll}, depth_id{1ll};

        for (auto id{0ll}; id < grid_z.mesh_size(); ++id)
        {
            CHECK(ring.density()(id) == density);
            CHECK(ring.specific_heat_capacity()(id) == capacity);
            CHECK(ring.heat_conductivity()(id) == conductivity);
            CHECK(ring.inner_radius(id) + ring.thickness(id) == ring.outer_radius(id));
            CHECK(ring.inner_radius(id) == 0.0);
            CHECK(ring.thickness(id) == ring.outer_radius(id));

            if (ring.thickness(id) == 0.0)
                CHECK(ring.outer_radius(id) == column.inner_radius(id));
            else
                CHECK(ring.outer_radius(id) == tube.inner_radius(id));
        }
    }
    {
        const auto &data2 = data["completion"]["variable"]["annulus"];

        const VarAnnulus annulus_ring_props{
            Completion::Density{data2["density"].get<VR>()},
            Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()},
            Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()}};

        casing.emplace(
            Completion::MaterialType::Annulus,
            FactoryVarRing::create_annulus_ring(
                annulus_ring_props,
                Completion::VarDepth{data2["depth_interval"].get<VR>()},
                casing.at(MaterialType::Tube),
                casing.at(MaterialType::Column),
                grid_z));
    }
    {
        const auto &data2 = data["completion"]["variable"]["cement"];
        casing.emplace(Completion::MaterialType::Cement,
                       VarRing{
                           VarRingSimple{
                               VarCement{
                                   Completion::Density{data2["density"].get<VR>()},
                                   Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<VR>()},
                                   Completion::HeatConductivity{data2["heat_conductivity"].get<VR>()}},
                               Completion::VarThickness{data2["radial_thickness"].get<VR>()},
                               Completion::VarDepth{data2["depth_interval"].get<VR>()}},
                           Completion::VarInnerRadius{data2["inner_radius"].get<VR>()}});
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
