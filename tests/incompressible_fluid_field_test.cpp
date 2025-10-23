
#include <vector>
#include <fstream>
#include <numbers>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Map/Grids2DMap.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/History/History.hpp>
#include <Injector/History/RatesFactory.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well/WellHoles.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Hydrodynamic/Incompressible/IncompressibleFluid.hpp>

#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/JT_FieldFactory.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/make_history.hpp"
#include "includes/set_is_permeable_stencils.hpp"
#include "includes/make_r_stencils.hpp"
#include "includes/get_completion.hpp"

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

const RealType tol = 1e-12;

TEST_CASE("Solver", "SelfSimilarCyl")
{
    ifstream f("heatflow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);

    /*fluid*/
    RealType
        viscosity{data["fluid"]["viscosity"]},
        density{data["fluid"]["density"]},
        capacity{data["fluid"]["specific_heat_capacity"]},
        heat_conductivity{data["fluid"]["heat_conductivity"]},
        joule_thomson{data["fluid"]["joule_thomson"]};

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
    CHECK(grid_rocks_z.mesh_size() == grid_z.mesh_size());
    CHECK(grid_rocks_r.mesh_size() == grid_r.mesh_size() - left_margin);
    for (auto row{0ll}; row < grid_rocks_r.mesh_size(); ++row)
    {
        CHECK(grid_rocks_r.mesh_nodes(row) == grid_r.mesh_nodes(row + left_margin));
        CHECK(grid_rocks_r.control_volumes(row) == grid_r.control_volumes(row + left_margin));
    }
    for (auto col{0ll}; col < grid_rocks_z.mesh_size(); ++col)
    {
        CHECK(grid_rocks_z.mesh_nodes(col) == grid_z.mesh_nodes(col));
        CHECK(grid_rocks_z.control_volumes(col) == grid_z.control_volumes(col));
    }

    Logs::Rocks::CoreSampleLogs
        core_logs{
            is_permeable_stencils,
            is_perforated_stencils,
            porosity_stencils,
            permeability_stencils,
            grid_z};

    Logs::Hydrodynamics::BaseHydrodynamics
        base_hydrodynamics(
            core_logs,
            medium_compressibility_stencils,
            ext_pressure_stencils);

    Properties::Rocks::RocksProps
        rock_field_props{
            base_hydrodynamics,
            grid2D_rocks};

    CHECK(rock_field_props.permeability_axes1.rows() == grid_z.mesh_size());
    CHECK(rock_field_props.permeability_axes1.cols() == grid_r.mesh_size() - left_margin);
    CHECK(rock_field_props.permeability_axes2.rows() == grid_z.mesh_size());
    CHECK(rock_field_props.permeability_axes2.cols() == grid_r.mesh_size() - left_margin);

    for (auto col{0ll}; col < grid_rocks_r.mesh_size(); ++col)
    {
        for (auto row{0ll}; row < grid_rocks_z.mesh_size(); ++row)
        {
            CHECK(rock_field_props.permeability_axes1.value(row, col) == 0.0);
            CHECK(rock_field_props.permeability_axes2.value(row, col) == core_logs.permeability(row));
        }
    }

    FaceProperties::Rocks::RocksFaceProps
        rock_face_props{
            rock_field_props,
            grid2D_rocks};
    CHECK(rock_face_props.permeability.face_vals_axes1.rows() == grid_rocks_z.dual_size() - 2ll);
    CHECK(rock_face_props.permeability.face_vals_axes1.cols() == grid_rocks_r.mesh_size());
    CHECK(rock_face_props.permeability.face_vals_axes2.rows() == grid_rocks_z.mesh_size());
    CHECK(rock_face_props.permeability.face_vals_axes2.cols() == grid_rocks_r.dual_size() - 2ll);

    for (auto col{0ll}; col < rock_face_props.permeability.face_vals_axes1.cols(); ++col)
    {
        for (auto row{0ll}; row < rock_face_props.permeability.face_vals_axes1.rows(); ++row)
        {
            CHECK(rock_face_props.permeability.face_vals_axes1(row, col) == 0.0);
        }
    }
    for (auto col{0ll}; col < rock_face_props.permeability.face_vals_axes2.cols(); ++col)
    {
        for (auto row{0ll}; row < rock_face_props.permeability.face_vals_axes1.rows(); ++row)
        {
            CHECK_THAT(rock_face_props.permeability.face_vals_axes2(row, col),
                       WithinRel(
                           core_logs.permeability(row) /
                               std::log(grid_rocks_r.mesh_nodes(col + 1ll) /
                                        grid_rocks_r.mesh_nodes(col)),
                           tol));
        }
    }

#pragma endregion

    const auto rMin{grid_r.dual_front()};
    const auto rMax{grid_r.dual_back()};

    // make fluid
    const PhasePropertiesJT water{
        FluidFactory::create_water_JT(
            Viscosity{viscosity},
            GPN::Density{density},
            GPN::SpecificHeatCapacity{capacity},
            GPN::HeatConductivity{heat_conductivity},
            JouleThomson{joule_thomson})};

    const auto RFP_weights{
        RFPFactory::create_from_container(
            RFP_weights_stencils,
            core_logs.is_permeable)};
    CHECK(RFP_weights.log_vals.sum() == 1.0);

    const CrossFlows cross_flows{
        RFP_weights, from_coords, to_layers};
    const auto WFP_weights{
        create_WFP(
            core_logs.is_perforated,
            RFP_weights,
            cross_flows)};
    CHECK(WFP_weights.log_vals.sum() == 1.0);

    const Well_Explicit well_explicit{
        core_logs.is_permeable, core_logs.is_perforated, RFP_weights};

    const Well_CrossFlow well{RFP_weights, WFP_weights, cross_flows};

    // history
    const shared_ptr<History> history{make_shared<History>(make_history(data))};

    using IncompressibleFluidField_t =
        decltype(IncompressibleFluidField{
            start_time,
            water,
            core_logs.permeability,
            base_hydrodynamics.ext_pressure,
            well,
            history,
            grid2D});

    auto ptr_pressure_field{
        make_shared<IncompressibleFluidField_t>(
            start_time,
            water,
            core_logs.permeability,
            base_hydrodynamics.ext_pressure,
            well,
            history,
            grid2D)};

    auto &pressure_field{*ptr_pressure_field};

    const RealType pi{std::numbers::pi_v<RealType>};
    const RealType tol{1e-10};

    for (auto i{0ll}; i < history->size(); ++i)
    {
        const auto r{history->get_record(i)};
        const auto rfp{well.get_RFP(r)};

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
                    CHECK(P.value(row, col) ==
                          base_hydrodynamics.ext_pressure(row) -
                              rfp(row) * water.viscosity /
                                  (2.0 * pi * core_logs.permeability(row) * grid_z.control_volumes(row)) *
                                  std::log(completion.sandface_radius(row) / rMax));
                }

                for (auto col{3ll}; col < grid_r.mesh_size(); ++col)
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
                for (auto col{3ll}; col < grid_r.mesh_size(); ++col)
                    CHECK(
                        P.value(row, col) == base_hydrodynamics.ext_pressure(row));
            }
        }
#pragma endregion
    }

    // rates field factory
    FaceProperties::IncompressibleRatesFactory rates_factory{
        ptr_pressure_field,
        grid2D, well, history, water};

    // mock SolverManaer behavior
    for (auto t{0ll}; t < history->size(); ++t)
    {
        history->advance();
        rates_factory.set_flow_field(history->time_moments[t], history->time_steps[t]);

        const auto &flux2_pos{rates_factory.get_heat_flow_in_axes2_pos()};
        const auto &flux2_neg{rates_factory.get_heat_flow_in_axes2_neg()};
        const auto &pressure{rates_factory.get_pressure_field().values()};

        const auto JT_term{Properties::JT_FieldFactory::create(rates_factory)};

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
        }
    }
}