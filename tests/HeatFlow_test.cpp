#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <numbers>
#include <cmath>
#include <vector>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Map/Grids2DMap.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/History/History.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Well/WellReservoirFlowProfileControl.hpp>
#include <Injector/Model/Well/WellBottomHolePressureControl.hpp>
#include <Injector/Model/Well/WellBottomHoleRateControl.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleRatesFactory.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluid.hpp>
#include <Injector/Model/Heat/HeatBoundaryConditions.hpp>
#include <Injector/Model/Completion.hpp>
#include <Injector/Model/ExtrudedCasingFactory.hpp>

#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/FlowField.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/make_r_stencils.hpp"
#include "includes/make_history.hpp"
#include "includes/make_geotherma.hpp"
#include "includes/make_water.hpp"
#include "includes/get_completion.hpp"
#include "includes/generate_stencils_and_steps.hpp"
#include "includes/set_is_permeable_stencils.hpp"
#include "includes/IC_BC.hpp"

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;
using json = nlohmann::json;

using namespace std;
using namespace GPN;
using namespace GPN::Logs;
using namespace GPN::CrossFlow;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::Completion;
using namespace GPN::Hydrodynamic;
using namespace GPN::Wells::ResFlowProfileControl;
using namespace GPN::Wells::BotHolePresControl;
using namespace GPN::Wells::BotHoleRateControl;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::FullImplicit;

TEST_CASE("Solver", "SelfSimilarCyl")
{
    ifstream f("heatflow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);

    /*START*/
    // input parameters
    /*fluid*/
    RealType
        viscosity{data["fluid"]["viscosity"].get<RealType>()},
        density{data["fluid"]["density"].get<RealType>()},
        capacity{data["fluid"]["specific_heat_capacity"].get<RealType>()},
        heat_conductivity{data["fluid"]["heat_conductivity"].get<RealType>()},
        joule_thomson{data["fluid"]["joule_thomson"].get<RealType>()};
    /*collector*/
    const auto thickness{data["collector"]["thickness"].get<VR>()};
    // const ptrdiff_t nLayers{thickness.size()};
    // hydrodynamic logs
    const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
    const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"].get<VR>())};
    const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>(), 1e-12)};
    const auto RFP_weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};
    const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
    const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};
    const auto is_permeable_stencils{set_is_permeable_stencils(is_perforated_stencils, to_layers)};
    const auto ext_pressure_stencils{transfer_to_eigen(data["collector"]["external_pressure"].get<VR>(), 1e5)};
    const auto medium_compressibility_stencils{transfer_to_eigen(data["collector"]["medium_compressibility"].get<VR>())};
    // heat logs
    const auto solid_heatconductivity_stencils{transfer_to_eigen(data["collector"]["heatConductivity"].get<VR>())};
    const auto solid_density_stencils{transfer_to_eigen(data["collector"]["solidDensity"].get<VR>())};
    const auto solid_specific_heatcapacity_stencils{transfer_to_eigen(data["collector"]["solidSpecificHeatCapacity"].get<VR>())};
    /*grid*/
    const auto
        z_minor_step{data["grid"]["z_minor_step"].get<RealType>()}; // m
    //  const ptrdiff_t rNodes{data["grid"]["rNodes"]};
    /*history*/
    const auto t_minor_step{read_minor_step(data)};
    const auto start_time{data["history"]["start_time"].get<RealType>()};
    /*temperatures*/
    /*completion*/
    // z-refiner
    const auto z_stencils{Grids::Factory::generate_dual_grid_stencils_from_steps(
        0.0, thickness)};
    RefinerVerticle z_refiner{z_minor_step, is_permeable_stencils};
    const auto temp_grid_z{
        Factory::create_axes<CoordinateTypes::Z>(
            z_refiner, z_stencils)};
    const Casing<VarRing> completion{get_completion(data, temp_grid_z)};
    /*END*/

    // make fluid
    const PhasePropertiesJT water{make_water(data)};

    // make grid2D
    // r_stencils
    const auto well_holes{WellHoles{WellHolesFactory::create(completion)}};
    const VR r_stencils{
        well_holes.get_stencils(
            data["grid"]["r_start"].get<RealType>(),
            data["grid"]["r_end"].get<RealType>())};
    // r-refiner
    const AbstractRefinerRadial *r_refiner{
        make_r_refiner(data)};

    const auto r_nodes{r_refiner->refine(r_stencils)};
    // the grid itself
    const auto grid2D{
        Grids::CylinderGridFactory::create(
            z_refiner, z_stencils,
            r_nodes)};
    const auto &grid_z{grid2D->first_coord()};
    const auto &grid_r{grid2D->second_coord()};

    const cptr<Grids::CylinderGridRock> grid2D_rocks{
        make_shared<Grids::CylinderGridRock>(grid2D)};
    constexpr auto left_margin{3ll};
    const auto &grid_rocks_z{grid2D_rocks->first_coord()};
    const auto &grid_rocks_r{grid2D_rocks->second_coord()};

    const auto rMin{grid_r.dual_front()};
    const auto rMax{grid_r.dual_back()};

    const ExtrudedCasing extr_completion{
        VarExtrudedCasingFactory::create(completion)};

    // cout << "radial dual grid stencils:\n"
    //      << grid_r.dual_nodes.transpose() << endl;

    // cout << "radial grid:\n"
    //      << grid2D->second_coord.dual_nodes.transpose() << endl;
    // cout << "vertical grid:\n"
    //      << grid2D->first_coord.dual_nodes.transpose() << endl;

    // cout << "radial grid cell centers:\n"
    //      << grid2D->second_coord.mesh_nodes.transpose() << endl;
    // cout << "vertical grid cell centers:\n"
    //      << grid2D->first_coord.mesh_nodes.transpose() << endl;

    // cout << "radial grid mesh steps:\n"
    //      << grid2D->second_coord.mesh_steps.transpose() << endl;
    // cout << "vertical grid mesh steps:\n"
    //      << grid2D->first_coord.mesh_steps.transpose() << endl;

    const Logs::Rocks::CoreSampleLogs core_logs{
        is_permeable_stencils,
        is_perforated_stencils,
        porosity_stencils,
        permeability_stencils,
        grid_z};

    const Logs::Hydrodynamics::BaseHydrodynamics base_hydrodynamics{
        Logs::Rocks::CoreSampleLogs{
            is_permeable_stencils,
            is_perforated_stencils,
            porosity_stencils,
            permeability_stencils,
            grid_z},
        medium_compressibility_stencils,
        ext_pressure_stencils};

    const Properties::Rocks::RocksProps
        rock_field_props{
            base_hydrodynamics,
            water,
            grid2D_rocks};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        solid_heatconductivity_stencils,
        porosity_stencils,
        water,
        grid_z};

    Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};
    heat_props.apply_well(extr_completion);

    FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};
    heat_face_props.apply_well(extr_completion);

    std::unique_ptr<const Logs::Geotherma> geotherma{
        make_unique<Logs::Geotherma>(
            make_geotherma(data, grid2D))};

    SECTION("WellReservoirFlowProfileControl")
    {
        const RFP_weights RFP_w{
            RFPFactory::create_from_container<Logs::RFP_weights>(
                StepProperty{RFP_weights_stencils},
                core_logs.is_permeable)};
        const CrossFlows cross_flows{
            from_coords, to_layers, RFP_w, core_logs.is_perforated};

        const shared_ptr<History> history{make_shared<History>(make_history(data))};
        using Well_t =
            decltype(WellReservoirFlowProfileControl{
                rock_field_props,
                cross_flows,
                history,
                grid2D_rocks});

        const ptr<Well_t> well{
            std::make_shared<Well_t>(
                rock_field_props,
                cross_flows,
                history,
                grid2D_rocks)};

        // fluid model for the pressure field
        using CompressibleFluidField_t =
            decltype(CompressibleFluidField{
                start_time,
                water,
                rock_field_props,
                *well,
                history,
                grid2D_rocks});

        auto ptr_pressure_field{
            make_shared<CompressibleFluidField_t>(
                start_time,
                water,
                rock_field_props,
                *well,
                history,
                grid2D_rocks)};

        // rates field factory
        using RatesFactory_t =
            decltype(FaceProperties::CompressibleRatesFactory{
                ptr_pressure_field,
                grid2D_rocks, well, history, water});
        auto ptr_rates_factory{make_shared<RatesFactory_t>(
            ptr_pressure_field,
            grid2D_rocks, well, history, water)};
        // initial condition
        const auto initial_state{GPN::ICFactory(start_time, grid2D, *geotherma)};
        // boundary conditions
        const GPN::Heat::HeatBC bc{
            grid2D,
            std::make_shared<GPN::FunctorBC<
                History,
                CompressibleFluidField_t,
                RatesFactory_t>>(
                history, ptr_rates_factory, *geotherma, grid2D),
            ptr_rates_factory};
        // solver

        using Solver_t = decltype(Solver{
            heat_face_props.medium_heat_conductivity,
            grid2D,
            heat_props.medium_vol_heatcapacity,
            ptr_rates_factory, initial_state,
            bc, start_time});

        auto solver_ptr{std::make_shared<Solver_t>(
            heat_face_props.medium_heat_conductivity,
            grid2D,
            heat_props.medium_vol_heatcapacity,
            ptr_rates_factory, initial_state,
            bc, start_time)};

        const auto &solver{*solver_ptr};

        SolverManager solver_manager{history, solver_ptr};

        solver_manager.run(t_minor_step);

        const auto &[p_times, p_states] = ptr_rates_factory->solution;

        // assert solution
        const double tol = 1E-11;
        const auto precision{1e-5};

        const auto &rates_factory{*ptr_rates_factory};
        // {
        //     string path{std::string{"flow_field.txt"}};
        //     ofstream f{path};
        //     f << (rates_factory.get_heat_flow_in_axes1() / precision).round() * precision << endl
        //       << endl;
        //     f << (rates_factory.get_heat_flow_in_axes2() / precision).round() * precision << endl
        //       << endl;
        //     f.close();
        // }
        // maximum principle
        const auto &[times, states] = solver.solution();
        // overall heat balance
        RealType cur_heat_incr = 0.0;
        RealType cum_inlet_heat = 0.0;
        //  cout << "volumetric heat capacity\n"
        //       << heat_props.medium_vol_heatcapacity.its_values << endl;
        for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
        {
            cur_heat_incr +=
                ((states[t].cur_state - states[t - 1ll].cur_state) *
                 heat_props.medium_vol_heatcapacity.values() * grid2D->volumes())
                    .sum();
            cum_inlet_heat +=
                (times[t] - times[t - 1ll]) *
                history->rates(t - 1ll) *
                water.volumetric_heat_capacity * (history->temps(t - 1ll) /*- initial_temperature*/);

            for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                cur_heat_incr +=
                    ((states[t].cur_state - states[t - 1ll].cur_state) *
                     heat_props.medium_vol_heatcapacity.values() * grid2D->volumes())
                        .sum();
                cum_inlet_heat +=
                    (times[t] - times[t - 1ll]) *
                    history->rates(t - 1ll) *
                    water.volumetric_heat_capacity * (history->temps(t - 1ll) /*- initial_temperature*/);

                RealType rel_tol = std::abs(2.0 * (cur_heat_incr - cum_inlet_heat) / (cur_heat_incr + cum_inlet_heat));
                //    CHECK(rel_tol < 0.05);
            }
        }
    }
    SECTION("WellBottomHolePressureControl")
    {
        data["history"]["control_type"] = "bottomhole_pressure";
        const CrossFlows cross_flows{
            from_coords, to_layers, core_logs.is_perforated};
        const shared_ptr<History> history{make_shared<History>(make_history(data))};
        using Well_t =
            decltype(WellBottomHolePressureControl{
                rock_field_props,
                cross_flows,
                history,
                water,
                grid2D_rocks});

        const ptr<Well_t> well{
            std::make_shared<Well_t>(
                rock_field_props,
                cross_flows,
                history,
                water,
                grid2D_rocks)};

        // fluid model for the pressure field
        using CompressibleFluidField_t =
            decltype(CompressibleFluidField{
                start_time,
                water,
                rock_field_props,
                *well,
                history,
                grid2D_rocks});

        auto ptr_pressure_field{
            make_shared<CompressibleFluidField_t>(
                start_time,
                water,
                rock_field_props,
                *well,
                history,
                grid2D_rocks)};

        // rates field factory
        using RatesFactory_t =
            decltype(FaceProperties::CompressibleRatesFactory{
                ptr_pressure_field,
                grid2D_rocks, well, history, water});
        auto ptr_rates_factory{make_shared<RatesFactory_t>(
            ptr_pressure_field,
            grid2D_rocks, well, history, water)};
        // initial condition
        const auto initial_state{GPN::ICFactory(start_time, grid2D, *geotherma)};
        // boundary conditions
        const GPN::Heat::HeatBC bc{
            grid2D,
            std::make_shared<GPN::FunctorBC<
                History,
                CompressibleFluidField_t,
                RatesFactory_t>>(
                history, ptr_rates_factory, *geotherma, grid2D),
            ptr_rates_factory};
        // solver

        using Solver_t = decltype(Solver{
            heat_face_props.medium_heat_conductivity,
            grid2D,
            heat_props.medium_vol_heatcapacity,
            ptr_rates_factory, initial_state,
            bc, start_time});

        auto solver_ptr{std::make_shared<Solver_t>(
            heat_face_props.medium_heat_conductivity,
            grid2D,
            heat_props.medium_vol_heatcapacity,
            ptr_rates_factory, initial_state,
            bc, start_time)};

        const auto &solver{*solver_ptr};

        SolverManager solver_manager{history, solver_ptr};

        solver_manager.run(t_minor_step);

        const auto &[p_times, p_states] = ptr_rates_factory->solution;

        // assert solution
        const double tol = 1E-11;
        const auto precision{1e-5};

        const auto &rates_factory{*ptr_rates_factory};
        // {
        //     string path{std::string{"flow_field.txt"}};
        //     ofstream f{path};
        //     f << (rates_factory.get_heat_flow_in_axes1() / precision).round() * precision << endl
        //       << endl;
        //     f << (rates_factory.get_heat_flow_in_axes2() / precision).round() * precision << endl
        //       << endl;
        //     f.close();
        // }
        // maximum principle
        const auto &[times, states] = solver.solution();
        // overall heat balance
        RealType cur_heat_incr = 0.0;
        RealType cum_inlet_heat = 0.0;
        //  cout << "volumetric heat capacity\n"
        //       << heat_props.medium_vol_heatcapacity.its_values << endl;
        for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
        {
            cur_heat_incr +=
                ((states[t].cur_state - states[t - 1ll].cur_state) *
                 heat_props.medium_vol_heatcapacity.values() * grid2D->volumes())
                    .sum();
            cum_inlet_heat +=
                (times[t] - times[t - 1ll]) *
                history->rates(t - 1ll) *
                water.volumetric_heat_capacity * (history->temps(t - 1ll) /*- initial_temperature*/);

            for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                cur_heat_incr +=
                    ((states[t].cur_state - states[t - 1ll].cur_state) *
                     heat_props.medium_vol_heatcapacity.values() * grid2D->volumes())
                        .sum();
                cum_inlet_heat +=
                    (times[t] - times[t - 1ll]) *
                    history->rates(t - 1ll) *
                    water.volumetric_heat_capacity * (history->temps(t - 1ll) /*- initial_temperature*/);

                RealType rel_tol = std::abs(2.0 * (cur_heat_incr - cum_inlet_heat) / (cur_heat_incr + cum_inlet_heat));
                //    CHECK(rel_tol < 0.05);
            }
        }
    }
    SECTION("WellBottomHoleRateControl")
    {
        data["history"]["control_type"] = "bottomhole_rate";
        const CrossFlows cross_flows{
            from_coords, to_layers, core_logs.is_perforated};
        const shared_ptr<History> history{make_shared<History>(make_history(data))};
        using Well_t =
            decltype(WellBottomHoleRateControl{
                rock_field_props,
                cross_flows,
                history,
                grid2D_rocks});

        const ptr<Well_t> well{
            std::make_shared<Well_t>(
                rock_field_props,
                cross_flows,
                history,
                grid2D_rocks)};

        // fluid model for the pressure field
        using FluidField_t =
            decltype(CompressibleFluidField{
                start_time,
                water,
                rock_field_props,
                *well,
                history,
                grid2D_rocks});

        auto ptr_pressure_field{
            make_shared<FluidField_t>(
                start_time,
                water,
                rock_field_props,
                *well,
                history,
                grid2D_rocks)};

        // rates field factory
        using RatesFactory_t =
            decltype(FaceProperties::CompressibleRatesFactory{
                ptr_pressure_field,
                grid2D_rocks, well, history, water});
        auto ptr_rates_factory{make_shared<RatesFactory_t>(
            ptr_pressure_field,
            grid2D_rocks, well, history, water)};
        // initial condition
        const auto initial_state{GPN::ICFactory(start_time, grid2D, *geotherma)};
        // boundary conditions
        const GPN::Heat::HeatBC bc{
            grid2D,
            std::make_shared<GPN::FunctorBC<
                History,
                FluidField_t,
                RatesFactory_t>>(
                history, ptr_rates_factory, *geotherma, grid2D),
            ptr_rates_factory};
        // solver

        using Solver_t = decltype(Solver{
            heat_face_props.medium_heat_conductivity,
            grid2D,
            heat_props.medium_vol_heatcapacity,
            ptr_rates_factory, initial_state,
            bc, start_time});

        auto solver_ptr{std::make_shared<Solver_t>(
            heat_face_props.medium_heat_conductivity,
            grid2D,
            heat_props.medium_vol_heatcapacity,
            ptr_rates_factory, initial_state,
            bc, start_time)};

        const auto &solver{*solver_ptr};

        SolverManager solver_manager{history, solver_ptr};

        solver_manager.run(t_minor_step);

        const auto &[p_times, p_states] = ptr_rates_factory->solution;

        // assert solution
        const double tol = 1E-11;
        const auto precision{1e-5};

        const auto &rates_factory{*ptr_rates_factory};
        // {
        //     string path{std::string{"flow_field.txt"}};
        //     ofstream f{path};
        //     f << (rates_factory.get_heat_flow_in_axes1() / precision).round() * precision << endl
        //       << endl;
        //     f << (rates_factory.get_heat_flow_in_axes2() / precision).round() * precision << endl
        //       << endl;
        //     f.close();
        // }
        // maximum principle
        const auto &[times, states] = solver.solution();
        // overall heat balance
        RealType cur_heat_incr = 0.0;
        RealType cum_inlet_heat = 0.0;
        //  cout << "volumetric heat capacity\n"
        //       << heat_props.medium_vol_heatcapacity.its_values << endl;
        for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
        {
            cur_heat_incr +=
                ((states[t].cur_state - states[t - 1ll].cur_state) *
                 heat_props.medium_vol_heatcapacity.values() * grid2D->volumes())
                    .sum();
            cum_inlet_heat +=
                (times[t] - times[t - 1ll]) *
                history->rates(t - 1ll) *
                water.volumetric_heat_capacity * (history->temps(t - 1ll) /*- initial_temperature*/);

            for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                cur_heat_incr +=
                    ((states[t].cur_state - states[t - 1ll].cur_state) *
                     heat_props.medium_vol_heatcapacity.values() * grid2D->volumes())
                        .sum();
                cum_inlet_heat +=
                    (times[t] - times[t - 1ll]) *
                    history->rates(t - 1ll) *
                    water.volumetric_heat_capacity * (history->temps(t - 1ll) /*- initial_temperature*/);

                RealType rel_tol = std::abs(2.0 * (cur_heat_incr - cum_inlet_heat) / (cur_heat_incr + cum_inlet_heat));
                //    CHECK(rel_tol < 0.05);
            }
        }
    }

}