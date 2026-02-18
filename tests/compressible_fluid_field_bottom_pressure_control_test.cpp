
#include <vector>
#include <fstream>
#include <numbers>
#include <cmath>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Map/Grids2DMap.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well/WellHoles.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Well/WellBottomHolePressureControl.hpp>

#include <Injector/Model/Hydrodynamic/Compressible/CompressibleRatesFactory.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluid.hpp>
#include <Injector/Model/Heat/HeatBoundaryConditions.hpp>

#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/JT_FieldFactory.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/make_history.hpp"
#include "includes/make_geotherma.hpp"
#include "includes/set_is_permeable_stencils.hpp"
#include "includes/make_r_stencils.hpp"
#include "includes/get_completion.hpp"
#include "includes/make_water.hpp"
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
using namespace GPN::Wells::BotHolePresControl;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::FullImplicit;

using VR = std::vector<GPN::RealType>;

const RealType pi{std::numbers::pi};
const RealType tol{1e-8};
const RealType exact_tol{1e-10};

TEST_CASE("Solver", "SelfSimilarCyl")
{
    ifstream f("heatflow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);
    data["history"]["control_type"] = "bottomhole_pressure";

    /*collector*/
    const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
    const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};

    const RealType zTop{0.0};
    const auto
        z_minor_step{data["grid"]["z_minor_step"].get<RealType>()}; // m
    const auto thickness{data["collector"]["thickness"].get<VR>()};
    const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>())};
    const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"].get<VR>())};
    const auto medium_compressibility_stencils{transfer_to_eigen(data["collector"]["medium_compressibility"].get<VR>())};

    const auto RFP_weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};
    const auto ext_pressure_stencils{transfer_to_eigen(data["collector"]["external_pressure"].get<VR>(), 1e5)};
    const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
    const auto is_permeable_stencils{set_is_permeable_stencils(is_perforated_stencils, to_layers)};

    // heat logs
    const auto solid_heatconductivity_stencils{transfer_to_eigen(data["collector"]["heatConductivity"].get<VR>())};
    const auto solid_density_stencils{transfer_to_eigen(data["collector"]["solidDensity"].get<VR>())};
    const auto solid_specific_heatcapacity_stencils{transfer_to_eigen(data["collector"]["solidSpecificHeatCapacity"].get<VR>())};

    const auto start_time{data["history"]["start_time"].get<RealType>()};

#pragma region MAKE-FLUID
    const PhasePropertiesJT water{make_water(data)};
#pragma endregion

    /*completion*/
    // z-refiner
    const auto z_stencils{Grids::Factory::generate_dual_grid_stencils_from_steps(
        0.0, thickness)};
    RefinerVerticle z_refiner{z_minor_step, is_permeable_stencils};
    const auto temp_grid_z{
        Factory::create_axes<CoordinateTypes::Z>(
            z_refiner, z_stencils)};
    const Casing<VarRing> completion{get_completion(data, temp_grid_z)};

    // make grid2D
    // r_stencils
    const auto well_holes{WellHoles{WellHolesFactory::create(completion)}};
    const auto r_well{well_holes.sandface_radius};
    const VR r_stencils{well_holes.get_stencils(
        data["grid"]["r_start"],
        data["grid"]["r_end"])};
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

    const ExtrudedCasing extr_completion{
        VarExtrudedCasingFactory::create(completion)};

    Logs::Rocks::CoreSampleLogs
        core_logs{
            is_permeable_stencils,
            is_perforated_stencils,
            porosity_stencils,
            permeability_stencils,
            grid_rocks_z};

    const auto &is_permeable{core_logs.is_permeable};
    const auto &permeability{core_logs.permeability};
    const auto &porosity{core_logs.porosity};
    const auto &cell_thickness{grid2D_rocks->first_coord().volumes()};

    Logs::Hydrodynamics::BaseHydrodynamics
        base_hydrodynamics(
            core_logs,
            medium_compressibility_stencils,
            ext_pressure_stencils);
    const auto &ext_pressure{base_hydrodynamics.ext_pressure};

    Properties::Rocks::RocksProps
        rock_field_props{
            base_hydrodynamics,
            water,
            grid2D_rocks};
    const auto &medium_compressibility_field{
        rock_field_props.medium_compressibility};
    CHECK(medium_compressibility_field.rows() == grid_rocks_z.mesh_size());
    CHECK(medium_compressibility_field.cols() == grid_rocks_r.mesh_size());

#pragma region HEAT-SETTINGS
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
#pragma endregion

#pragma region MAKE-HISTORY
    const ptr<History> history{make_shared<History>(make_history(data))};
#pragma endregion
#pragma region MAKE-WELL
    const CrossFlows cross_flows{
        from_coords, to_layers, core_logs.is_perforated};

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

    const auto &PI{well->PI};
    for (auto row{0ll}; row < grid_rocks_z.mesh_size(); ++row)
    {
        const auto r3{grid_rocks_r.mesh_nodes(0ll)};
        const auto sandface{grid_rocks_r.dual_nodes(0ll)};
        const auto val{
            2.0 * pi * cell_thickness(row) * permeability(row) / water.viscosity /
            std::log(r3 / sandface)};
        CHECK_THAT(PI(row), WithinRel(val, exact_tol));
        CHECK(PI(row) >= 0.0);
    }

#pragma endregion
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
    auto &pressure_field{*ptr_pressure_field};

    // rates field factory
    using RatesFactory_t =
        decltype(FaceProperties::CompressibleRatesFactory{
            ptr_pressure_field,
            grid2D_rocks, well, history, water});
    auto ptr_rates_factory{make_shared<RatesFactory_t>(
        ptr_pressure_field,
        grid2D_rocks, well, history, water)};
    RatesFactory_t &rates_factory{*ptr_rates_factory};

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

    //    const auto rMin{grid_r.dual_front()};
    const auto rMax{grid_r.dual_back()};

    const auto mu{water.viscosity};
    const auto r_sandface{well_holes.sandface_radius};
    CHECK(r_sandface == grid_r.dual_nodes(3ll));
    CHECK(r_sandface == grid_rocks_r.dual_front());
    const auto &time_intervals{history->time_steps};
    const auto numerical_step{read_minor_step(data)};
    RealType cur_time{start_time};
    ptrdiff_t counter{0ll};
    // mock SolverManager::run
    for (auto t_step{0ll}; t_step < history->time_steps.size(); ++t_step)
    {
        history->advance();

        const size_t internal_step_count{
            static_cast<size_t>(
                std::abs(std::ceil(time_intervals[t_step] / numerical_step)))};
        const RealType step{time_intervals[t_step] / internal_step_count};

        for (size_t id{0ull}; id < internal_step_count; ++id, cur_time += step)
        {
            const auto t{cur_time + step};
            const auto collector_pressure_prev{ptr_pressure_field->get_rock_pressure()};
            rates_factory.set_flow_field(cur_time, step);
            const auto record{history->get_current_record()};
            const auto collector_pressure{ptr_pressure_field->get_rock_pressure()};
            const RealType P_bot{well->P_bot(record, collector_pressure)};
            //        std::cout << "P_bot:        " << P_bot << std::endl;
            //        std::cout << "P_bot_solver: " << P_bot_solver << std::endl;
            CHECK_THAT(P_bot, WithinRel(record.pressure, exact_tol));
#pragma region VERIFY-PRESSURE-PROBLEM-MATRIX
            {
                const auto &A{pressure_field.get_solver()->get_problem_matrix()};

                const auto nz{grid_rocks_z.mesh_size()};
                const auto nr{grid_rocks_r.mesh_size()};
                for (auto row{0ll}; row < grid_rocks_z.mesh_size(); ++row)
                {
                    bool flag{(is_permeable(row) == 1.0) || (is_permeable(row) == 0.0)};
                    REQUIRE(flag);
#pragma region VERIFY-PERMEABLE-LAYERS
                    if (is_permeable(row) == 1.0)
                    {
#pragma region SANDFACE-BOUNDARY
                        {
                            const auto col{0ll};
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            for (auto idx{0ll}; idx < l; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                            {
                                const auto idx{l};
                                const auto val{
                                    grid2D_rocks->volume(row, col) *
                                        medium_compressibility_field.value(row, col) / step +
                                    2.0 * pi * permeability(row) * cell_thickness(row) / water.viscosity *
                                        (1.0 / std::log(grid_rocks_r.mesh_nodes(col + 1ll) / grid_rocks_r.mesh_nodes(col))) +
                                    PI(row)};

                                CHECK_THAT(A.coeff(l, idx),
                                           WithinRel(val, exact_tol));
                            }

                            for (auto idx{l + 1ll}; idx < l + nz; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                            {
                                const auto idx{l + nz};
                                const auto val{
                                    -2.0 * pi * permeability(row) * cell_thickness(row) / water.viscosity *
                                    (1.0 / std::log(grid_rocks_r.mesh_nodes(col + 1ll) / grid_rocks_r.mesh_nodes(col)))};

                                CHECK_THAT(A.coeff(l, idx),
                                           WithinRel(val, exact_tol));
                            }

                            for (auto idx{l + nz + 1ll}; idx < A.cols(); ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                        }
#pragma endregion
                        for (auto col{1ll}; col < grid_rocks_r.mesh_size() - 1ll; ++col)
                        {
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            for (auto idx{0ll}; idx < l - nz; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                            {
                                const auto idx{l - nz};
                                const auto val{
                                    -2.0 * pi * permeability(row) * cell_thickness(row) / water.viscosity *
                                    (1.0 / std::log(grid_rocks_r.mesh_nodes(col) / grid_rocks_r.mesh_nodes(col - 1ll)))};

                                CHECK_THAT(A.coeff(l, idx),
                                           WithinRel(val, exact_tol));
                            }
                            for (auto idx{l - nz + 1ll}; idx < l; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);

                            const auto val{
                                grid2D_rocks->volume(row, col) *
                                    medium_compressibility_field.value(row, col) / step +
                                2.0 * pi * permeability(row) * cell_thickness(row) / water.viscosity *
                                    (1.0 / std::log(grid_rocks_r.mesh_nodes(col + 1ll) / grid_rocks_r.mesh_nodes(col)) +
                                     1.0 / std::log(grid_rocks_r.mesh_nodes(col) / grid_rocks_r.mesh_nodes(col - 1ll)))};

                            CHECK_THAT(A.coeff(l, l),
                                       WithinRel(val, exact_tol));

                            for (auto idx{l + 1ll}; idx < l + nz; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                            {
                                const auto idx{l + nz};
                                const auto val{
                                    -2.0 * pi * permeability(row) * cell_thickness(row) / water.viscosity *
                                    (1.0 / std::log(grid_rocks_r.mesh_nodes(col + 1ll) / grid_rocks_r.mesh_nodes(col)))};
                                CHECK_THAT(A.coeff(l, idx),
                                           WithinRel(val, exact_tol));
                            }
                            for (auto idx{l + nz + 1ll}; idx < A.cols(); ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                        }
#pragma region EXTERNAL-DOMAIN-BOUNDARY
                        {
                            const auto col{grid_rocks_r.mesh_size() - 1ll};
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            for (auto idx{0ll}; idx < l; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                            CHECK(A.coeff(l, l) == 1.0);
                            for (auto idx{l + 1ll}; idx < A.cols(); ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                        }
#pragma endregion
                    }
#pragma endregion
#pragma region VERIFY-NON_PERMEABLE-LAYERS
                    else
                    {
                        for (auto col{0ll}; col < grid_rocks_r.mesh_size() - 1ll; ++col)
                        {
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            for (auto idx{0ll}; idx < l; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                            CHECK_THAT(A.coeff(l, l),
                                       WithinRel(grid2D_rocks->volume(row, col) / step, exact_tol));
                            for (auto idx{l + 1ll}; idx < A.cols(); ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                        }
#pragma region EXTERNAL-DOMAIN-BOUNDARY
                        {
                            const auto col{grid_rocks_r.mesh_size() - 1ll};
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            for (auto idx{0ll}; idx < l; ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                            CHECK(A.coeff(l, l) == 1.0);
                            for (auto idx{l + 1ll}; idx < A.cols(); ++idx)
                                CHECK(A.coeff(l, idx) == 0.0);
                        }
#pragma endregion
                    }
#pragma endregion
                }
            }
#pragma endregion

#pragma region VERIFY-PRESSURE-PROBLEM-RHS
            {
                const auto &rhs{pressure_field.get_solver()->get_problem_rhs()};
                for (auto row{0ll}; row < grid_rocks_z.mesh_size(); ++row)
                {
                    if (is_permeable(row) == 1.0)
                    {
                        { // sandface boundary
                            const auto col{0ll};
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            const auto val{collector_pressure_prev(row, col) * grid2D_rocks->volume(row, col) *
                                               medium_compressibility_field.value(row, col) / step +
                                           P_bot * PI(row)};
                            CHECK_THAT(rhs(l),
                                       WithinRel(val, exact_tol));
                        }
                        for (auto col{1ll}; col < grid_rocks_r.mesh_size() - 1ll; ++col)
                        {
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            const auto val{collector_pressure_prev(row, col) * grid2D_rocks->volume(row, col) *
                                           medium_compressibility_field.value(row, col) / step};
                            CHECK_THAT(rhs(l),
                                       WithinRel(val, exact_tol));
                        }
                        { // external contour
                            const auto col{grid_rocks_r.mesh_size() - 1ll};
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            CHECK(rhs(l) == ext_pressure(row));
                        }
                    }
                    else
                    {
                        for (auto col{0ll}; col < grid_rocks_r.mesh_size(); ++col)
                        {
                            const auto l{grid2D_rocks->to_linear(row, col)};
                            CHECK(rhs(l) == 0.0);
                            CHECK(rhs(l) == ext_pressure(row));
                        }
                    }
                }
            }
#pragma endregion
            const auto rfp{well->get_RFP(history->get_current_record())};
            const auto Q{rfp.sum()};
#pragma region CHECK-PRESSURE
            const auto &P{pressure_field.current_pressure().values()};
            CHECK(rfp.rows() == grid_z.mesh_size());

            for (auto row{0ll}; row < grid_z.mesh_size(); ++row)
            {
                const auto rate{rfp(row)};
                const auto p_ex{ext_pressure(row)};
                const auto h{grid_z.control_volumes(row)};
                const auto k{permeability(row)};
                const auto beta{base_hydrodynamics.medium_compressibility(row)};
#pragma region CHECK-RFP
                CHECK(rate == -PI(row) * (collector_pressure(row, 0ll) - P_bot));
#pragma endregion

                // pressure is const inside completion
                const auto P_sandface{well->get_pressure_at_sandface(record, collector_pressure)};
                CHECK(P_sandface(row) == P(row, 2ll));
                for (auto col{0ll}; col < 2ll; ++col)
                {
                    CHECK(P(row, col) == P(row, 2ll));
                }
                if (is_permeable(row) == 1.0)
                {     // in permeable layers
                    { // in the cement
                        auto col{2ll};
                        CHECK(grid_r.dual_nodes(col + 1ll) == completion.sandface_radius(row));
                        INFO("col: " << col);
                        CHECK_THAT(P(row, col),
                                   WithinRel(
                                       P(row, col + 1ll) + rfp(row) * mu /
                                                               (2.0 * pi * k * h) *
                                                               std::log(grid_r.mesh_nodes(3ll) / completion.sandface_radius(row)),
                                       tol));
                    }
                    // external contour
                    {
                        const auto col{grid_r.mesh_size() - 1ll};
                        const auto r{grid_r.mesh_nodes(col)};
                        const auto calc_val{P(row, col) - p_ex};
                        INFO("row: " << row << "; r: " << r);
                        CHECK_THAT(
                            calc_val / 1e5,
                            WithinAbs(0.0,
                                      exact_tol * 100));
                    }
                }
                else
                { // in impermeable layers
                    for (auto col{0ll}; col < grid_r.mesh_size(); ++col)
                        CHECK(
                            P(row, col) == ext_pressure(row));
                }
#pragma region CECK-WELL-PRESSURE

#pragma endregion
            }
#pragma endregion
#pragma region CHECK-RATES
            const auto C{water.volumetric_heat_capacity};
            const auto flux1{rates_factory.get_heat_flow_in_axes1_neg() + rates_factory.get_heat_flow_in_axes1_pos()};
            CHECK(flux1.rows() == grid_z.dual_size());
            CHECK(flux1.cols() == grid_r.mesh_size());
            const auto flux2{rates_factory.get_heat_flow_in_axes2_neg() + rates_factory.get_heat_flow_in_axes2_pos()};
            CHECK(flux2.rows() == grid_z.mesh_size());
            CHECK(flux2.cols() == grid_r.dual_size());
            const auto &mobility2{pressure_field.face_mobility.face_vals_axes2};
            const auto wfp{well->get_WFP(history->get_current_record())};
            //        std::cout << "wfp:\n" << wfp.transpose() << std::endl;
            CHECK(wfp.rows() == grid_z.mesh_size());
            CHECK_THAT(rfp.sum(), WithinRel(wfp.sum(), exact_tol));
            const auto cement_flow{well->get_verticle_cement_flow(history->get_current_record())};
            CHECK(cement_flow.rows() == grid_z.dual_size());
            const auto well_flow{well->get_verticle_well_flow(history->get_current_record())};
            CHECK(well_flow.rows() == grid_z.dual_size());
            RealType well_loss_cum_sum{0.0};
            RealType well_accum_cum_sum{0.0};
#pragma region HORIZONTAL-RATES
            for (auto row{0ll}; row < grid_z.mesh_size(); ++row)
            {
                const auto k{permeability(row)};
                const auto h{grid_z.volume(row)};

                if (is_permeable(row) == 1.0)
                { // check horizontal rates in permeable layer
                    {
                        const auto col{0ll};
                        CHECK(flux2(row, col) == 0.0);
                    }
                    {
                        for (auto col{1ll}; col <= 2ll; ++col)
                        {
                            CHECK(flux2(row, col) == C * wfp(row));
                        }
                    }
                    {
                        const auto col{left_margin};
                        const auto r{grid_r.mesh_nodes(col)};
                        const auto ref_rate{
                            -2.0 * pi * k * h / mu *
                            (P(row, col) - P(row, col - 1ll)) / std::log(r / r_sandface)};
                        const auto rate{rfp(row)};
                        CHECK_THAT(ref_rate,
                                   WithinRel(
                                       rate /*mobility2(row, col - left_margin)*/,
                                       1e-4));
                        CHECK_THAT(flux2(row, col),
                                   WithinRel(C * rate, tol));
                    }
                    for (auto col{left_margin + 1ll}; col < grid_z.dual_size() - 1ll; ++col)
                    {
                        const auto r{grid_r.mesh_nodes(col - 1ll)};
                        const auto r_next{grid_r.mesh_nodes(col)};
                        const auto ref_rate{
                            -2.0 * C * pi * k * h / mu *
                            (P(row, col) - P(row, col - 1ll)) / std::log(r_next / r)};
                        CHECK_THAT(flux2(row, col),
                                   WithinRel(ref_rate, 1e-10));
                    }
                    {
                        const auto col{grid_r.dual_size() - 1ll};
                        CHECK(flux2(row, col) == 0.0);
                    }
                    // check rate signs
                    for (auto col{left_margin + 1ll}; col < grid_z.dual_size() - 1ll; ++col)
                    {
                        const auto r{grid_r.mesh_nodes(col - 1ll)};
                        const auto r_next{grid_r.mesh_nodes(col)};
                        const auto ref_rate{
                            -2.0 * C * pi * k * h / mu *
                            (P(row, col) - P(row, col - 1ll)) / std::log(r_next / r)};
                        CHECK(flux2(row, col) * flux2(row, left_margin) >= 0.0);
                    }
                }
                else
                { // check horizontal rates in impermeable layers
                    {
                        const auto col{0ll};
                        CHECK(flux2(row, col) == 0.0);
                    }
                    for (auto col{1ll}; col < left_margin; ++col)
                    {
                        // this comes from cross flow data
                        INFO("row: " << row << ", col: " << col);
                        CHECK(flux2(row, col) == C * wfp(row));
                    }

                    for (auto col{left_margin}; col < grid_r.dual_size(); ++col)
                    {
                        INFO("row: " << row << ", col: " << col);
                        CHECK(flux2(row, col) == 0.0);
                    }
                }

                if (core_logs.is_permeable(row) == 1.0)
                {
                    if (base_hydrodynamics.medium_compressibility(row) == 0.0)
                    {
                        const auto r_max{grid_r.dual_back()};
                        const auto thickness{grid2D_rocks->first_coord().volume(row)};
                        const auto k{permeability(row)};
                        const auto p_ext{ext_pressure(row)};
                        for (auto col{3ll}; col < grid_r.mesh_size() - 1ll; ++col)
                        {
                            const auto r{grid_r.mesh_nodes(col)};
                            const RealType ref_p{p_ext -
                                                 rfp(row) * water.viscosity / (2.0 * pi * thickness * k) * std::log(r / r_max)};
                            CHECK_THAT(ref_p, WithinRel(P(row, col), tol));
                        }
                    }
                }
            }
#pragma endregion
#pragma region VERTICAL-RATES
            //    std::cout << "well_loss_cum_sum:\n" << well_loss_cum_sum.transpose() << std::endl;
            //    std::cout << "wfp:\n" << wfp.transpose() << std::endl;
            //    std::cout << "Q: " << Q << std::endl;

            for (auto row{grid_z.dual_size() - 1ll}; row >= 0ll; --row)
            {
                { // flow in the tube
                    const auto col{0ll};
                    INFO("row: " << row);
                    if (well_accum_cum_sum == 0.0)
                    {
                        CHECK_THAT(
                            well_accum_cum_sum,
                            WithinAbs(well_flow(row), exact_tol));
                    }
                    else
                    {
                        CHECK_THAT(well_accum_cum_sum,
                                   WithinRel(well_flow(row), tol));
                    }

                    if (std::abs(flux1(row, col) / C / Q) < exact_tol)
                    {
                        INFO("C: " << C << ", acc_flux: " << well_accum_cum_sum << ", flux1: " << flux1(row, col));
                        CHECK_THAT(well_accum_cum_sum,
                                   WithinAbs(flux1(row, col) / C, exact_tol));
                    }
                    else
                    {
                        INFO("C: " << C << ", acc_flux: " << well_accum_cum_sum << ", flux1: " << flux1(row, col));
                        CHECK_THAT(well_accum_cum_sum,
                                   WithinRel(flux1(row, col) / C, tol));
                    }

                    if (row > 0ll)
                        well_accum_cum_sum += wfp(row - 1ll);
                }
                { // flow in the sandwich
                    const auto col{1ll};
                    INFO("row: " << row);
                    CHECK(flux1(row, col) == 0.0);
                }
                { // flow in the cement
                    const auto col{2ll};
                    INFO("row: " << row);
                    CHECK(flux1(row, col) / C == cement_flow(row));
                }
                for (auto col{3ll}; col < grid_r.mesh_size(); ++col)
                {
                    INFO("row: " << row);
                    CHECK(flux1(row, col) == 0.0);
                }
            }
#pragma endregion
#pragma endregion
#pragma region CHECK-JT
            const auto &flux2_pos{rates_factory.get_heat_flow_in_axes2_pos()};
            const auto &flux2_neg{rates_factory.get_heat_flow_in_axes2_neg()};
            const auto &pressure{P};

            const auto JT_term{Properties::JT_FieldFactory::create_spatial(rates_factory)};

            const auto JT_temporal_term{Properties::JT_FieldFactory::create_temporal(rates_factory)};

            CHECK(JT_term.rows() == JT_temporal_term.rows());
            CHECK(JT_term.cols() == JT_temporal_term.cols());

            for (auto row{0ll}; row < JT_term.rows(); ++row)
            {
                {
                    const auto col{0ll};
                    CHECK_THAT(
                        water.JT * (flux2_pos(row, col) * (pressure(row, col)) +
                                    flux2_neg(row, col + 1ll) * (pressure(row, col + 1ll) - pressure(row, col))),
                        WithinRel(JT_term.value(row, col), tol));
                }
                for (auto col{1ll}; col < JT_term.cols() - 1ll; ++col)
                {
                    // INFO("row: " << row << "col: " << col);
                    // CHECK_THAT(
                    //     water.JT*(flux2_pos(row, col) * (pressure(row, col) - pressure(row, col - 1ll)) +
                    //         flux2_neg(row, col + 1ll) * (pressure(row, col + 1ll) - pressure(row, col))),
                    //     WithinRel(JT_term.value(row, col), tol));

                    INFO("row: " << row << "col: " << col);
                    CHECK_THAT(
                        water.JT * ((flux2_pos(row, col) - flux2_neg(row, col + 1ll)) * pressure(row, col) +
                                    flux2_neg(row, col + 1ll) * pressure(row, col + 1ll) -
                                    flux2_pos(row, col) * pressure(row, col - 1ll)),
                        WithinRel(JT_term.value(row, col), tol));
                }
                {
                    const auto col{JT_term.cols() - 1ll};
                    // CHECK_THAT(
                    //     water.JT*(flux2_pos(row, col) * (pressure(row, col) - pressure(row, col - 1ll)) +
                    //         flux2_neg(row, col + 1ll) * (-pressure(row, col))),
                    //     WithinRel(JT_term.value(row, col), tol));

                    CHECK_THAT(
                        water.JT * ((flux2_pos(row, col) - flux2_neg(row, col + 1ll)) * pressure(row, col) -
                                    flux2_pos(row, col) * pressure(row, col - 1ll)),
                        WithinRel(JT_term.value(row, col), tol));
                }

                {
                    for (auto col{0ll}; col < left_margin; ++col)
                        CHECK(JT_temporal_term.value(row, col) == 0.0);
                }
                {
                    for (auto col{left_margin}; col < JT_temporal_term.cols(); ++col)
                    {
                        const auto ref{porosity(row) * water.adiabatic_factor * grid2D->volume(row, col) *
                                       (rates_factory.pressure_field->P->value(row, col) -
                                        rates_factory.pressure_field->P_prev->value(row, col))};
                        INFO("row: " << row << ", col: " << col);
                        CHECK(ref == JT_temporal_term.value(row, col));
                    }
                }
            }

#pragma endregion
        }
    }
    //    INFO("At least a single point should be checked by analytical expression!");
    //    CHECK(counter > 0ll);
}