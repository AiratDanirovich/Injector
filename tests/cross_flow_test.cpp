#include <fstream>
#include <iostream>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/History/History.hpp>
#include <Injector/History/RatesFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Hydrodynamic/Incompressible/IncompressibleFluid.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/set_is_permeable_stencils.hpp"
#include "includes/make_history.hpp"

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;
using VR = std::vector<GPN::RealType>;

using namespace std;
using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::Hydrodynamic;
using namespace GPN::FaceProperties;
using namespace GPN::Logs;
using namespace GPN::CrossFlow;

RealType viscosity{6e-4}, density{1000}, capacity{4200}, heat_conductivity{0.6};

TEST_CASE("CrossFlow", "")
{
    ifstream f("heatflow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);
    // hydrodynamic logs
    const auto RFP_weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};
    const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
    const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};
    const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
    const auto is_permeable_stencils{set_is_permeable_stencils(is_perforated_stencils, to_layers)};
    const auto ext_pressure_stencils{transfer_to_eigen(data["collector"]["external_pressure"].get<VR>(), 1e5)};
    const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>(), 1e-12)};
    /*collector*/
    const VR thickness{data["collector"]["thickness"].get<VR>()};
    const auto start_time{data["history"]["start_time"].get<RealType>()};
    /*grid*/
    const RealType
        z_minor_step{data["grid"]["z_minor_step"]}; // m

    const auto grid_z{
        Factory::create_axes<CoordinateTypes::Z>(
            RefinerVerticle{
                data["grid"]["z_minor_step"].get<RealType>(),
                is_permeable_stencils},
            Grids::Factory::generate_dual_grid_stencils_from_steps(
                0.0, data["collector"]["thickness"].get<VR>()))};

    const auto is_permeable{
        IsPermeableFactory::create(is_permeable_stencils, grid_z)};
    const auto is_perforated{
        IsPermeableFactory::create(is_perforated_stencils, grid_z)};
    const auto permeability{
        Logs::PermeabilityFactory::create(
            permeability_stencils,
            is_permeable_stencils,
            grid_z)};

    const auto RFP_weights{
        RFPFactory::create_from_container(
            RFP_weights_stencils,
            is_permeable)};

    const RealType total_rate{1.0};
    const CrossFlows cross_flows{
        RFP_weights, from_coords, to_layers};

    cout << "flux weights: " << RFP_weights.log_vals.transpose() << endl;
    const auto &dual_nodes{RFP_weights.grid.dual_nodes};
    const auto &dual_stencils{RFP_weights.grid.dual_stencils};
    for (auto i{0ull}; i < cross_flows.cross_flow_data.size(); ++i)
    {
        const auto &cf{cross_flows.cross_flow_data[i]};
        cout << "flux vector: " << cf.verticle_flux.transpose() << endl;

        auto to_id{0ll};
        { // get TO id
            const auto to_coord{dual_stencils(to_layers[i] + 1ll)};
            while (dual_nodes(to_id) < to_coord)
                ++to_id;
            --to_id;
            CHECK(dual_nodes(to_id) < to_coord);
            CHECK(dual_nodes(to_id + 1ll) == to_coord);
        }
        auto from_id{0ll};
        { // get FROM id
            while (dual_nodes(from_id) < from_coords[i])
                ++from_id;
            --from_id;
            CHECK(dual_nodes(from_id) < from_coords[i]);
            CHECK(dual_nodes(from_id + 1ll) > from_coords[i]);
        }

        const RealType dir{from_id < to_id ? 1.0 : -1.0};
        cout << "dir: " << dir << std::endl;
        for (auto j{0ll}; j < std::min(to_id, from_id) + 1ll; ++j)
        {
            INFO("top_id: " << j << ", to_id: " << to_id << ", from_id: " << from_id);
            CHECK(cf.verticle_flux(j) == 0.0);
        }
        for (auto j{std::max(to_id, from_id) + 1ll}; j < dual_nodes.rows(); ++j)
            CHECK(cf.verticle_flux(j) == 0.0);

        for (auto j{std::min(to_id, from_id) + 1ll}; j < std::max(to_id, from_id) + 1ll; ++j)
        {
            INFO("top_id: " << j << ", to_id: " << to_id << ", from_id: " << from_id << ", dir: " << dir);
            CHECK(cf.verticle_flux(j) == dir * RFP_weights(to_id));
        }
    }

    const auto WFP_weights{
        create_WFP(
            is_perforated,
            RFP_weights,
            cross_flows)};

    const Well_CrossFlow well{RFP_weights, WFP_weights, cross_flows};

    // make fluid
    const PhaseProperties water{
        FluidFactory::create_water(
            Viscosity{viscosity},
            GPN::Density{density},
            GPN::SpecificHeatCapacity{capacity},
            GPN::HeatConductivity{heat_conductivity})};

    // history
    const History history{make_history(data)};

    // z-refiner
    RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
    const VR r_stencils{0.0, 0.1, 0.2, 0.3};
    const auto grid2D{
        Grids::CylinderGridFactory::create(refiner,
                                           Grids::Factory::generate_dual_grid_stencils_from_steps(
                                               0.0, thickness),
                                           r_stencils)};
    // external pressure log
    const auto external_pressure{
        Logs::ExtPressureFactory::create(
            ext_pressure_stencils,
            is_permeable_stencils,
            grid_z)};
    // fluid model for the pressure field
    auto pressure_field{
        IncompressibleFluidField(
            start_time,
            water,
            permeability,
            external_pressure,
            well,
            grid2D)};
    // rates field factory
    FaceProperties::IncompressibleRatesFactory rates_factory{
        pressure_field,
        grid2D, well, history, water};

    rates_factory.set_flow_field(0.0, 1800.0);
}
