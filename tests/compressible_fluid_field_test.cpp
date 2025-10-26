
#include <vector>
#include <fstream>
#include <numbers>
#include <cmath>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Map/Grids2DMap.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/History/RatesFactory.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well/WellHoles.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluid.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluidSolver.hpp>

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

using VR = std::vector<GPN::RealType>;

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
    const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>(), 1e-12)};
    const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"].get<VR>())};
    const auto medium_compressibility_stencils{transfer_to_eigen(data["collector"]["medium_compressibility"].get<VR>())};

    const auto RFP_weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};
    const auto ext_pressure_stencils{transfer_to_eigen(data["collector"]["external_pressure"].get<VR>(), 1e5)};
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

    const cptr<Grids::CylinderGridRock> grid2D_rocks{
        make_shared<Grids::CylinderGridRock>(grid2D)};
    constexpr auto left_margin{3ll};
    const auto &grid_rocks_z{grid2D_rocks->first_coord()};
    const auto &grid_rocks_r{grid2D_rocks->second_coord()};

    Logs::Rocks::CoreSampleLogs
        core_logs{
            is_permeable_stencils,
            is_perforated_stencils,
            porosity_stencils,
            permeability_stencils,
            grid_rocks_z};

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
#pragma region MAKE-WELL
    const auto RFP_weights{
        RFPFactory::create_from_container(
            RFP_weights_stencils,
            core_logs.is_permeable)};
    const CrossFlows cross_flows{
        RFP_weights, from_coords, to_layers};
    const auto WFP_weights{
        create_WFP(
            core_logs.is_perforated,
            RFP_weights,
            cross_flows)};
    const Well_CrossFlow well{RFP_weights, WFP_weights, cross_flows};
#pragma endregion
#pragma region MAKE-HISTORY
    const ptr<History> history{make_shared<History>(make_history(data))};
#pragma endregion
    using CompressibleFluidField_t =
        decltype(CompressibleFluidField{
            start_time,
            water,
            core_logs.permeability,
            rock_field_props,
            base_hydrodynamics.ext_pressure,
            well,
            history,
            grid2D_rocks});

    auto ptr_pressure_field{
        make_shared<CompressibleFluidField_t>(
            start_time,
            water,
            core_logs.permeability,
            rock_field_props,
            base_hydrodynamics.ext_pressure,
            well,
            history,
            grid2D_rocks)};

    auto &pressure_field{*ptr_pressure_field};

    const RealType pi{std::numbers::pi_v<RealType>};
    const RealType tol{1e-10};

    const auto &grid_z{grid2D->first_coord()};
    const auto &grid_r{grid2D->second_coord()};

    //    const auto rMin{grid_r.dual_front()};
    const auto rMax{grid_r.dual_back()};

    const auto mu{water.viscosity};
    const auto &time_intervals{history->time_steps};
    const auto numerical_step{data["history"]["t_minor_step"].get<RealType>()};
    RealType cur_time{start_time};
    for (auto t_step{0ll}; t_step < history->time_steps.size(); ++t_step)
    {
        history->advance();

        const auto record{history->get_current_record()};
        const auto rfp{well.get_RFP(record)};

        const size_t internal_step_count{
            static_cast<size_t>(
                std::abs(std::ceil(time_intervals[t_step] / numerical_step)))};
        const RealType step{time_intervals[t_step] / internal_step_count};

        for (size_t id{0ull}; id < internal_step_count; ++id)
        {
            cur_time += step;
            const auto t{cur_time};

            pressure_field.set_pressure_field(step, record);

            const auto &P{pressure_field.current_pressure()};

            for (auto row{0ll}; row < grid_z.mesh_size(); ++row)
            {
                // for (auto col{0ll}; col < 2ll; ++col)
                //     CHECK(P.value(row, col) == P.value(row, 2ll));
                if (core_logs.is_permeable(row) == 1.0)
                {
                    // {
                    //     auto col{2ll};
                    //     CHECK(grid_r.dual_nodes(col + 1ll) == completion.sandface_radius(row));
                    //     CHECK(P.value(row, col) ==
                    //           base_hydrodynamics.ext_pressure(row) -
                    //               rfp(row) * water.viscosity /
                    //                   (2.0 * pi * core_logs.permeability(row) * grid_z.control_volumes(row)) *
                    //                   std::log(completion.sandface_radius(row) / rMax));
                    // }

                    const auto h{grid_z.control_volumes(row)};
                    const auto k{core_logs.permeability(row)};
                    const auto beta{base_hydrodynamics.medium_compressibility(row)};
                    const auto piezo_cond{k / (mu * beta)};
                    if (r_well * r_well < 0.01 * piezo_cond * t)
                    {
                        const auto rate{rfp(row)};
                        const auto p_ex{base_hydrodynamics.ext_pressure(row)};

                        for (auto col{left_margin}; col < grid_r.mesh_size() - 1ll; ++col)
                        {
                            const auto r{grid_r.mesh_nodes(col)};
                            const auto ei{-std::expint(-r * r / (4.0 * piezo_cond * t))};
                            const auto ref_val{
                                //    p_ex -
                                rate * mu /
                                (4.0 * piezo_cond * k * h) * ei};
                            const auto calc_val{P.value(row, col)};

                            const auto stat_p{-rate * mu /
                                              (2.0 * pi * k * h) *
                                              std::log(r / rMax)};

                            INFO("row: " << row << "; col: " << col << "; r: " << r << "; ratio: " << stat_p / (calc_val - p_ex));
                            CHECK_THAT(
                                (calc_val - p_ex),
                                WithinRel(stat_p,
                                          tol));
                        }
                        {
                            const auto col{grid_r.mesh_size() - 1ll};
                            const auto r{grid_r.mesh_nodes(col)};
                            const auto calc_val{P.value(row, col)};
                            INFO("row: " << row << "; r: " << r);
                            CHECK_THAT(
                                (calc_val - p_ex) / 1e5,
                                WithinAbs(0.0,
                                          tol));
                        }
                    }
                }
                else
                {
                    for (auto col{left_margin}; col < grid_r.mesh_size(); ++col)
                        CHECK(
                            P.value(row, col) == base_hydrodynamics.ext_pressure(row));
                }
            }
        }
    }

    // rates field factory
    // FaceProperties::IncompressibleRatesFactory rates_factory{
    //     ptr_pressure_field,
    //     grid2D_rocks, well, history, water};

    // for (auto t{0ll}; t < history->size(); ++t)
    // { // mock SolverManager::run
    //     history->advance();
    //     rates_factory.set_flow_field(history->time_moments[t], history->time_steps[t]);

    //     const auto &flux2_pos{rates_factory.get_heat_flow_in_axes2_pos()};
    //     const auto &flux2_neg{rates_factory.get_heat_flow_in_axes2_neg()};
    //     const auto &pressure{rates_factory.get_pressure_field().values()};

    //     const auto JT_term{Properties::JT_FieldFactory::create(rates_factory)};

    //     for (auto row{0ll}; row < JT_term.rows(); ++row)
    //     {
    //         {
    //             const auto col{0ll};
    //             CHECK_THAT(
    //                 flux2_pos(row, col) * (pressure(row, col)) +
    //                     flux2_neg(row, col + 1ll) * (pressure(row, col + 1ll) - pressure(row, col)),
    //                 WithinRel(JT_term.value(row, col) / water.JT, tol));
    //         }
    //         for (auto col{1ll}; col < JT_term.cols() - 1ll; ++col)
    //         {
    //             CHECK_THAT(
    //                 flux2_pos(row, col) * (pressure(row, col) - pressure(row, col - 1ll)) +
    //                     flux2_neg(row, col + 1ll) * (pressure(row, col + 1ll) - pressure(row, col)),
    //                 WithinRel(JT_term.value(row, col) / water.JT, tol));
    //         }
    //         {
    //             const auto col{JT_term.cols() - 1ll};
    //             CHECK_THAT(
    //                 flux2_pos(row, col) * (pressure(row, col) - pressure(row, col - 1ll)) +
    //                     flux2_neg(row, col + 1ll) * (-pressure(row, col)),
    //                 WithinRel(JT_term.value(row, col) / water.JT, tol));
    //         }
    //     }
    // }
}