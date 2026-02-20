
#include <vector>
#include <fstream>
#include <numbers>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Map/Grids2DMap.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/History/History.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well/WellHoles.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Well/WellReservoirFlowProfileControl.hpp>
#include <Injector/Model/Hydrodynamic/Incompressible/IncompressibleRatesFactory.hpp>
#include <Injector/Model/Hydrodynamic/Incompressible/IncompressibleFluid.hpp>

#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/JT_FieldFactory.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/make_history.hpp"
#include "includes/set_is_permeable_stencils.hpp"
#include "includes/make_r_stencils.hpp"
#include "includes/get_completion.hpp"
#include "includes/make_water.hpp"

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

const RealType tol{1e-12};

TEST_CASE("Solver", "SelfSimilarCyl")
{
    ifstream f("heatflow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);

    /*collector*/
    const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
    const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};

    const RealType zTop{0.0};
    const auto
        z_minor_step{data["grid"]["z_minor_step"].get<RealType>()}; // m
    const auto thickness{data["collector"]["thickness"].get<VR>()};
    const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"].get<VR>())};
    const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>(), 1e-12)};
    const auto RFP_weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};
    const auto ext_pressure_stencils{transfer_to_eigen(data["collector"]["external_pressure"].get<VR>(), 1e5)};
    const auto medium_compressibility_stencils{transfer_to_eigen(data["collector"]["medium_compressibility"].get<VR>())};
    const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
    const auto is_permeable_stencils{set_is_permeable_stencils(is_perforated_stencils, to_layers)};

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
    const VR r_stencils{
        WellHoles{WellHolesFactory::create(completion)}.get_stencils(
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

#pragma region CHECK-MAP-GRID
    constexpr std::ptrdiff_t left_margin{3ll};
    const cptr<Grids::CylinderGridRock> grid2D_rocks{
        make_shared<Grids::CylinderGridRock>(grid2D)};
    const auto &grid_rocks_z{grid2D_rocks->first_coord()};
    const auto &grid_rocks_r{grid2D_rocks->second_coord()};

    for (auto col{0ll}; col < grid_rocks_r.mesh_size(); ++col)
    {
        for (auto row{0ll}; row < grid_rocks_z.mesh_size(); ++row)
        {
            CHECK(grid2D_rocks->volume(row, col) == grid2D->volume(row, col + left_margin));
            CHECK(grid2D_rocks->to_linear(row, col) == grid2D->to_linear(row, col));
        }
    }

    for (auto col{0ll}; col < grid_rocks_r.mesh_size(); ++col)
    {
        CHECK(
            grid2D_rocks->face_area_axes1(col) ==
            grid2D->face_area_axes1(col + left_margin));
    }
    for (auto row{0ll}; row < grid_rocks_z.mesh_size(); ++row)
    {
        CHECK(
            grid2D_rocks->face_area_axes2(row) ==
            grid2D->face_area_axes2(row));
    }

    CHECK(grid_rocks_z.mesh_size() == grid_z.mesh_size());
    CHECK(grid_rocks_z.dual_front() == grid_z.dual_front());
    CHECK(grid_rocks_z.dual_back() == grid_z.dual_back());
    CHECK(grid_rocks_z.dual_size() == grid_z.dual_size());

    CHECK(grid_rocks_r.mesh_size() == grid_r.mesh_size() - left_margin);
    CHECK(grid_rocks_r.dual_front() == grid_r.dual_nodes(left_margin));
    CHECK(grid_rocks_r.dual_back() == grid_r.dual_back());
    CHECK(grid_rocks_r.dual_size() == grid_r.dual_size() - left_margin);

    for (auto row{0ll}; row < grid_rocks_r.mesh_size(); ++row)
    {
        CHECK(grid_rocks_r.mesh_nodes(row) == grid_r.mesh_nodes(row + left_margin));
        CHECK(grid_rocks_r.dual_nodes(row) == grid_r.dual_nodes(row + left_margin));
        CHECK(grid_rocks_r.control_volumes(row) == grid_r.control_volumes(row + left_margin));
    }
    for (auto col{0ll}; col < grid_rocks_z.mesh_size(); ++col)
    {
        CHECK(grid_rocks_z.mesh_nodes(col) == grid_z.mesh_nodes(col));
        CHECK(grid_rocks_z.dual_nodes(col) == grid_z.dual_nodes(col));
        CHECK(grid_rocks_z.control_volumes(col) == grid_z.control_volumes(col));
    }

    Logs::Rocks::CoreSampleLogs
        core_logs{
            is_permeable_stencils,
            is_perforated_stencils,
            porosity_stencils,
            permeability_stencils,
            grid_z};
    const auto& porosity{core_logs.porosity};

    Logs::Hydrodynamics::BaseHydrodynamics
        base_hydrodynamics(
            core_logs,
            medium_compressibility_stencils,
            ext_pressure_stencils);

    Properties::Rocks::RocksProps
        rock_field_props{
            base_hydrodynamics,
            water,
            grid2D_rocks};

    CHECK(rock_field_props.mobility_axes1.rows() == grid_z.mesh_size());
    CHECK(rock_field_props.mobility_axes1.cols() == grid_r.mesh_size() - left_margin);
    CHECK(rock_field_props.mobility_axes2.rows() == grid_z.mesh_size());
    CHECK(rock_field_props.mobility_axes2.cols() == grid_r.mesh_size() - left_margin);

    for (auto col{0ll}; col < grid_rocks_r.mesh_size(); ++col)
    {
        for (auto row{0ll}; row < grid_rocks_z.mesh_size(); ++row)
        {
            CHECK(rock_field_props.mobility_axes1.value(row, col) == 0.0);
            CHECK_THAT(rock_field_props.mobility_axes2.value(row, col),
                       WithinRel(core_logs.permeability(row) / water.viscosity, tol));
        }
    }

    FaceProperties::Rocks::RocksFaceProps
        rock_face_props{
            rock_field_props,
            grid2D_rocks};
    CHECK(rock_face_props.mobility.face_vals_axes1.rows() == grid_rocks_z.dual_size() - 2ll);
    CHECK(rock_face_props.mobility.face_vals_axes1.cols() == grid_rocks_r.mesh_size());
    CHECK(rock_face_props.mobility.face_vals_axes2.rows() == grid_rocks_z.mesh_size());
    CHECK(rock_face_props.mobility.face_vals_axes2.cols() == grid_rocks_r.dual_size() - 2ll);

    for (auto col{0ll}; col < rock_face_props.mobility.face_vals_axes1.cols(); ++col)
    {
        for (auto row{0ll}; row < rock_face_props.mobility.face_vals_axes1.rows(); ++row)
        {
            CHECK(rock_face_props.mobility.face_vals_axes1(row, col) == 0.0);
        }
    }
    for (auto col{0ll}; col < rock_face_props.mobility.face_vals_axes2.cols(); ++col)
    {
        for (auto row{0ll}; row < rock_face_props.mobility.face_vals_axes1.rows(); ++row)
        {
            CHECK_THAT(rock_face_props.mobility.face_vals_axes2(row, col),
                       WithinRel(
                           core_logs.permeability(row) / water.viscosity /
                               std::log(grid_rocks_r.mesh_nodes(col + 1ll) /
                                        grid_rocks_r.mesh_nodes(col)),
                           tol));
        }
    }

#pragma endregion

    const auto rMin{grid_r.dual_front()};
    const auto rMax{grid_r.dual_back()};

    const Logs::RFP_weights RFP_w{
        RFPFactory::create_from_container<Logs::RFP_weights>(
            StepProperty{RFP_weights_stencils},
            core_logs.is_permeable)};
    CHECK(RFP_w.log_vals.sum() == 1.0);

    const CrossFlows cross_flows{
        from_coords, to_layers, RFP_w, core_logs.is_perforated};
    const auto WFP_weights{
        create_WFP_weights(
            core_logs.is_perforated,
            RFP_w,
            cross_flows.cross_flow_handler)};
    CHECK(WFP_weights.log_vals.sum() == 1.0);


    // history
    const shared_ptr<History> history{make_shared<History>(make_history(data))};

    using Well_t =
        decltype(WellReservoirFlowProfileControl{
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

    using FluidField_t =
        decltype(IncompressibleFluidField{
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

    auto &pressure_field{*ptr_pressure_field};

    const RealType pi{std::numbers::pi_v<RealType>};
    const RealType tol{1e-10};

    for (auto i{0ll}; i < history->size(); ++i)
    {
        const auto r{history->get_record(i)};
        const auto rfp{well->get_RFP(r)};

        pressure_field.set_pressure_field(0.0, r);

#pragma region VERIFY-PRESSURE
        const auto &P{pressure_field.current_pressure()};
        for (auto row{0ll}; row < grid_z.mesh_size(); ++row)
        {
            for (auto col{0ll}; col < 2ll; ++col)
                CHECK(P.value(row, col) == P.value(row, 2ll));
            if (core_logs.is_permeable(row) == 1.0)
            {
                {
                    auto col{2ll};
                    CHECK(grid_r.dual_nodes(col + 1ll) == completion.sandface_radius(row));
                    CHECK_THAT(
                        P.value(row, col),
                        WithinRel(
                            base_hydrodynamics.ext_pressure(row) -
                                rfp(row) * water.viscosity /
                                    (2.0 * pi * core_logs.permeability(row) * grid_z.control_volumes(row)) *
                                    std::log(completion.sandface_radius(row) / rMax),
                            tol));
                }

                for (auto col{left_margin}; col < grid_r.mesh_size(); ++col)
                    CHECK_THAT(
                        P.value(row, col),
                        WithinRel(
                            base_hydrodynamics.ext_pressure(row) -
                                rfp(row) * water.viscosity /
                                    (2.0 * pi * core_logs.permeability(row) * grid_z.control_volumes(row)) *
                                    std::log(grid_r.mesh_nodes(col) / rMax),
                            tol));
            }
            else
            {
                for (auto col{left_margin}; col < grid_r.mesh_size(); ++col)
                    CHECK(
                        P.value(row, col) == base_hydrodynamics.ext_pressure(row));
            }
        }
#pragma endregion
#pragma region VERIFY-RATES
        {
#pragma region CHECKS

            // TODO: fix volumetric balance of fluid flow in reservoir, and in the well (if possible)

            //     const auto rate{history->rate()};
            //     const auto pressure{history->pressure()};
            //     // verify flow field
            //     const auto v1{rates_factory.get_heat_flow_in_axes1()}; // verticle flow
            //     for (auto row{0ll}; row < v1.rows(); ++row)
            //     {
            //         CHECK(v1(row, 0ll) * rate >= 0.0); // flow in tube
            //         CHECK(v1(row, 1ll) == 0.0);        // flow in tube-wall::annulus::column-wall
            //         // CHECK(col == 2ll) // flow in cement is pretty much complex
            //         for (auto col{3ll}; col < v1.cols(); ++col)
            //             CHECK(v1(row, col) == 0.0);
            //     }
            //     const auto v2{rates_factory.get_heat_flow_in_axes2()}; // horizontal flow
            //     // boundary conditions at r = 0.0
            //     for (auto row{0ll}, col{0ll}; row < v2.rows(); ++row)
            //     {
            //         CHECK(v2(row, col) == 0.0);
            //         CHECK_THAT(v1(row, col), WithinRel(v1(row + 1, col) + v2(row, col + 1ll), tol));
            //     }
            //     // compare against WFP
            //     const auto WFP{(water.volumetric_heat_capacity * well.get_WFP(history->get_current_record())).eval()};
            //     for (auto row{0ll}, col{1ll}; row < v2.rows(); ++row)
            //     {
            //         CHECK(v2(row, col) == v2(row, 2ll));
            //         CHECK(v2(row, col) == WFP(row));
            //     }
            //     // compare against RFP
            //     const auto RFP{(water.volumetric_heat_capacity * well.get_RFP(history->get_current_record())).eval()};
            //     for (auto row{0ll}; row < v2.rows(); ++row)
            //     {
            //         CHECK(v2(row, 3ll) == RFP(row));
            //         for (auto col{4ll}; col < v2.cols(); ++col)
            //             CHECK(v2(row, col) == v2(row, 3ll));
            //     }

            //     // flow volume balance
            //     for (auto row{0ll}; row < grid2D->first_coord().mesh_size(); ++row)
            //     {
            //         for (auto col{0ll}; col < grid2D->second_coord().mesh_size(); ++col)
            //         {
            //             INFO("" << "col: " << col << ", row: " << row << ", bottom: " << -v1(row + 1ll, col) << ", top: " << v1(row, col) << ", right: " << v2(row, col + 1ll) << ", left: " << -v2(row, col));
            //             CHECK_THAT(-v1(row + 1ll, col) + v1(row, col), WithinRel(v2(row, col + 1ll) - v2(row, col), tol));
            //         }
            //     }
            // }
#pragma endregion
        }
#pragma endregion
    }

    // rates field factory
    FaceProperties::IncompressibleRatesFactory rates_factory{
        ptr_pressure_field,
        grid2D_rocks, well, history, water};

    // mock SolverManaer behavior
    for (auto t{0ll}; t < history->size(); ++t)
    {
        history->advance();
        rates_factory.set_flow_field(history->time_moments[t], history->time_steps[t]);

        const auto &flux2_pos{rates_factory.get_heat_flow_in_axes2_pos()};
        const auto &flux2_neg{rates_factory.get_heat_flow_in_axes2_neg()};
        const auto &pressure{rates_factory.get_pressure_field().values()};

        const auto JT_term{Properties::JT_FieldFactory::create_spatial(rates_factory)};
        const auto JT_temporal_term{Properties::JT_FieldFactory::create_temporal(rates_factory)};

        CHECK(JT_term.rows() == JT_temporal_term.rows());
        CHECK(JT_term.cols() == JT_temporal_term.cols());

        for (auto row{0ll}; row < JT_term.rows(); ++row)
        {
            if (core_logs.permeability(row) == 0.0)
                for (auto col{0ll}; col < JT_term.cols(); ++col)
                    CHECK(JT_term.value(row, col) == 0.0);
            {
                const auto col{0ll};
                CHECK_THAT(
                    water.JT * (flux2_pos(row, col) * (pressure(row, col)) +
                                flux2_neg(row, col + 1ll) * (pressure(row, col + 1ll) - pressure(row, col))),
                    WithinRel(JT_term.value(row, col), tol));
            }
            for (auto col{1ll}; col < JT_term.cols() - 1ll; ++col)
            {
                CHECK_THAT(
                    water.JT * (flux2_pos(row, col) * (pressure(row, col) - pressure(row, col - 1ll)) +
                                flux2_neg(row, col + 1ll) * (pressure(row, col + 1ll) - pressure(row, col))),
                    WithinRel(JT_term.value(row, col), tol));
            }
            {
                const auto col{JT_term.cols() - 1ll};
                CHECK_THAT(
                    water.JT * (flux2_pos(row, col) * (pressure(row, col) - pressure(row, col - 1ll)) +
                                flux2_neg(row, col + 1ll) * (-pressure(row, col))),
                    WithinRel(JT_term.value(row, col), tol));
            }

            {
                for (auto col{0ll}; col < left_margin; ++col)
                    CHECK(JT_temporal_term.value(row, col) == 0.0);
            }
            {
                for (auto col{left_margin}; col < JT_temporal_term.cols(); ++col)
                {
                    const auto ref{water.adiabatic_factor*grid2D->volume(row, col)*porosity(row) *
                                   (rates_factory.pressure_field->P->value(row, col) -
                                    rates_factory.pressure_field->P_prev->value(row, col))};
                    INFO(
                        "row: " << row << 
                        ", col: " << col << 
                        ", t_id: " << t << 
                        ", rows: " << grid2D_rocks->first_coord().mesh_size() << 
                        ", cols: " << grid2D->second_coord().mesh_size());
                    CHECK(ref == JT_temporal_term.value(row, col));
                }
            }
        }
    }
}