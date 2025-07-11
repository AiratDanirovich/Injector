#include <iostream>
#include <fstream>
#include <vector>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Model/WellFactory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include "includes/get_completion.hpp"
#include "includes/transfer_to_eigen.hpp"
#include "includes/make_r_stencils.hpp"
#include "includes/get_dist.hpp"

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;

using VR = std::vector<GPN::RealType>;

using namespace std;
using namespace Catch;
using namespace Catch::Matchers;
using namespace GPN;
using namespace GPN::Logs;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::Completion;

const RealType tol = 1e-12;

TEST_CASE("apply_well_test", "SelfSimilarCyl")
{
  ifstream f("apply_well_test_data.json");
  REQUIRE(f.is_open());
  json data = json::parse(f);

  /*START*/
  // input parameters
  /*fluid*/
  RealType
      viscosity{data["fluid"]["viscosity"]},
      density{data["fluid"]["density"]},
      capacity{data["fluid"]["specific_heat_capacity"]},
      heat_conductivity{data["fluid"]["heat_conductivity"]};
  /*collector*/
  const VR thickness = data["collector"]["thickness"];
  // hydrodynamic logs
  const auto is_permeable_stencils{transfer_to_eigen(data["collector"]["is_permeable"])};
  const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"])};
  const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"])};
  const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"], 1e-12)};
  const auto weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"])};
  // heat logs
  const auto solid_heatconductivity_stencils{transfer_to_eigen(data["collector"]["heatConductivity"])};
  const auto solid_density_stencils{transfer_to_eigen(data["collector"]["solidDensity"])};
  const auto solid_specific_heatcapacity_stencils{transfer_to_eigen(data["collector"]["solidSpecificHeatCapacity"])};
  /*grid*/
  const RealType
      z_minor_step{data["grid"]["z_minor_step"]}; // m
  /*completion*/
  const Casing completion{get_completion(data)};
  /*END*/

  // make grid2D
  // r_stencils
  const WellHoles well_holes{WellHolesFactory::create(completion)};
  const VR r_stencils{make_r_stencils(data, well_holes)};
  const RealType &rMax = r_stencils.back();
  const RealType &rMin = r_stencils.front();

  // CHECK r_stencils
  {
    CHECK(r_stencils[0ull] == 0.0);
    CHECK(r_stencils[1ull] == completion.flow_radius);
    CHECK(r_stencils[2ull] == completion.column_outer_radius);
    CHECK(r_stencils[3ull] == completion.sandface_radius);
  }

  // z-refiner
  RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
  // the grid itself
  const auto grid2D{
      Grids::CylinderGridFactory::create(refiner,
                                         Grids::Factory::generate_dual_grid_stencils_from_steps(
                                             0.0, thickness),
                                         r_stencils)};
  const auto &grid_z{grid2D->first_coord};
  const auto &grid_r{grid2D->second_coord};

  // cout << "radial dual grid stencils:\n"
  //      << transfer_to_eigen(r_stencils).transpose() << endl;

  const Logs::Rocks::CoreSampleLogs core_data{
      is_permeable_stencils,
      is_perforated_stencils,
      porosity_stencils,
      permeability_stencils,
      grid_z};

  // CHECK core_data
  {
    const auto &is_perm = core_data.is_permeable.log_vals;
    const auto &is_perf = core_data.is_perforated.log_vals;
    const auto &perm = core_data.permeability.log_vals;
    const auto &poro = core_data.porosity.log_vals;
    for (auto row{0ll}; row < is_perm.rows(); ++row)
    {
      // find val in coarse grid
      const ptrdiff_t idx{
          get_dist(
              grid_z.dual_stencils.dual_nodes,
              grid_z.mesh_nodes(row))};
      // solid density
      CHECK_THAT(is_perm(row), WithinRel(is_permeable_stencils(idx), tol));
      // solid_specific_heatcapacity
      CHECK_THAT(is_perf(row), WithinRel(is_perforated_stencils(idx), tol));
      // solid_heat_conductivity
      CHECK_THAT(perm(row), WithinRel(permeability_stencils(idx), tol));
      // solid_heat_conductivity
      CHECK_THAT(poro(row), WithinRel(porosity_stencils(idx), tol));
    }
  }

  // make fluid
  const PhaseProperties water{
      FluidFactory::create_water(
          Viscosity{viscosity},
          Density{density},
          SpecificHeatCapacity{capacity},
          GPN::HeatConductivity{heat_conductivity})};

  const auto weights{
      RateWeightsFactory::create(
          weights_stencils,
          core_data.is_permeable,
          grid_z)};

  // CHECK weights
  {
    const auto &w = weights.log_vals;
    for (auto row{0ll}; row < w.rows(); ++row)
    {
      // find val in coarse grid
      const ptrdiff_t idx{
          get_dist(
              grid_z.dual_stencils.dual_nodes,
              grid_z.mesh_nodes(row))};
      // solid density
      CHECK_THAT(w(row), WithinRel(weights_stencils(idx), tol));
    }
  }

  const Well_Explicit well{
      core_data.is_permeable, core_data.is_perforated, weights};

  const Logs::Rocks::HeatLogs heat_logs{
      solid_density_stencils,
      solid_specific_heatcapacity_stencils,
      solid_heatconductivity_stencils,
      porosity_stencils,
      water,
      grid2D->first_coord};

  // CHECK heat_logs
  {
    const auto &capacity = heat_logs.medium_vol_heatcapacity.log_vals;
    const auto &heat_conductivity = heat_logs.medium_heat_conductivity.log_vals;
    const auto &s_vol_heatcapacity = heat_logs.solid_vol_heatcapacity.log_vals;

    const auto &s_heat_conductivity = heat_logs.solid_heat_conductivity.log_vals;
    const auto &s_density = heat_logs.solid_density.log_vals;
    const auto &s_spesific_heat_cap = heat_logs.solid_specific_heatcapacity.log_vals;
    const auto &porosity = core_data.porosity.log_vals;
    for (auto row{0ll}; row < capacity.rows(); ++row)
    {
      // medium_vol_heatcapacity
      CHECK_THAT(capacity(row),
                 WithinRel(
                     porosity(row) * water.volumetric_heat_capacity +
                         (1 - porosity(row)) * s_vol_heatcapacity(row),
                     tol));
      // medium_heat_conductivity
      CHECK_THAT(heat_conductivity(row),
                 WithinRel(
                     porosity(row) * water.heat_conductivity +
                         (1 - porosity(row)) * s_heat_conductivity(row),
                     tol));
      // solid_vol_heatcapacity
      CHECK_THAT(s_vol_heatcapacity(row),
                 WithinRel(
                     s_density(row) * s_spesific_heat_cap(row),
                     tol));
      // find val in coarse grid
      const ptrdiff_t idx{
          get_dist(
              grid_z.dual_stencils.dual_nodes,
              grid_z.mesh_nodes(row))};
      // solid density
      CHECK_THAT(s_density(row), WithinRel(solid_density_stencils(idx), tol));
      // solid_specific_heatcapacity
      CHECK_THAT(s_spesific_heat_cap(row), WithinRel(solid_specific_heatcapacity_stencils(idx), tol));
      // solid_heat_conductivity
      CHECK_THAT(s_heat_conductivity(row), WithinRel(solid_heatconductivity_stencils(idx), tol));
    }
  }

  Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  // CHECK heat_props
  {
    const auto &capacity = heat_props.medium_vol_heatcapacity.values();
    const auto &heat_conductivity_1 = heat_props.medium_heat_conductivity_axes1.values();
    const auto &heat_conductivity_2 = heat_props.medium_heat_conductivity_axes2.values();
    for (auto row{0ll}; row < capacity.rows(); ++row)
    {
      for (auto col{0ll}; col < capacity.cols(); ++col)
      {
        CHECK_THAT(capacity(row, col),
                   WithinRel(
                       heat_logs.medium_vol_heatcapacity(row), tol));
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       heat_logs.medium_heat_conductivity(row), tol));
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       heat_conductivity_2(row, col), tol));
      }
    }
  }

  // properties of material that fills the well up to the Sandface
  heat_props.apply_well(completion, well);

  const auto &Tube{completion[MaterialType::Tube]};
  const auto &Annulus{completion[MaterialType::Annulus]};
  const auto &Column{completion[MaterialType::Column]};
  const auto &Sandface{completion[MaterialType::Cement]};

  const auto casing_vert_cond{std::numbers::pi * (Tube.heat_conductivity * Tube.thickness * (Tube.inner_radius + Tube.outer_radius) + Annulus.heat_conductivity * Annulus.thickness * (Annulus.inner_radius + Annulus.outer_radius) + Column.heat_conductivity * Column.thickness * (Column.inner_radius + Column.outer_radius)) /
                              (std::numbers::pi * (Column.outer_radius + Tube.inner_radius) * (Column.outer_radius - Tube.inner_radius))};
  CHECK_THAT(casing_vert_cond, WithinRel(completion.integral_vertical_casing_heat_conductivity(), tol));
  
  const auto cement_vert_cond{std::numbers::pi * (
    Sandface.heat_conductivity * Sandface.thickness * (Sandface.inner_radius + Sandface.outer_radius)) /
                              (std::numbers::pi * (Sandface.outer_radius + Sandface.inner_radius) * (Sandface.thickness))};
  CHECK_THAT(cement_vert_cond, WithinRel(completion.integral_vertical_cement_heat_conductivity(), tol));

  // CHECK heat_props --- after "apply_well"
  {
    const auto &capacity = heat_props.medium_vol_heatcapacity.values();
    const auto &heat_conductivity_1 = heat_props.medium_heat_conductivity_axes1.values();
    const auto &heat_conductivity_2 = heat_props.medium_heat_conductivity_axes2.values();
    for (auto row{0ll}; row < capacity.rows(); ++row)
    {
      { // col == 0
        const ptrdiff_t col = 0ll;
        // flow heat capacity is equal to fluid-water heat capacity
        CHECK_THAT(capacity(row, col),
                   WithinRel(
                       water.volumetric_heat_capacity, tol));
        CHECK_THAT(capacity(row, col),
                   WithinRel(
                       completion.front().volumetric_heat_capacity, tol));
        // vertical heat conductivity is equal to water
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       water.heat_conductivity, tol));
        // radial heat conductivity of flowing water is infinity
        CHECK(std::isinf(heat_conductivity_2(row, col)));
      }
      { // col == 1
        const ptrdiff_t col = 1ll;
        CHECK_THAT(capacity(row, col),
                   WithinRel(
                       (Tube.linear_heat_capacity +
                        Annulus.linear_heat_capacity +
                        Column.linear_heat_capacity +
                        Sandface.linear_heat_capacity) /
                           completion.area(),
                       tol));
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       casing_vert_cond
                       //  (Tube.integral_vertical_heat_conductivity +
                       //   Annulus.integral_vertical_heat_conductivity +
                       //   Column.integral_vertical_heat_conductivity +
                       //   CementInner.integral_vertical_heat_conductivity +
                       //   Sandface.integral_vertical_heat_conductivity) /
                       //      completion.area()
                       ,
                       tol));
        CHECK_THAT(heat_conductivity_2(row, col),
                   WithinRel(
                       Sandface.heat_conductivity /
                           std::log(Sandface.outer_radius / Sandface.inner_radius) *
                           std::log(grid_r.dual_nodes(2ll) / grid_r.mesh_nodes(1ll)),
                       tol));
      }

      for (auto col{2ll}; col < capacity.cols(); ++col)
      {
        CHECK_THAT(capacity(row, col),
                   WithinRel(
                       heat_logs.medium_vol_heatcapacity(row), tol));
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       heat_logs.medium_heat_conductivity(row), tol));
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       heat_conductivity_2(row, col), tol));
      }
    }
  }

  FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};

  // CHECK heat_face_props
  {
    const auto &f_conductivity_1 = heat_face_props.medium_heat_conductivity.face_vals_axes1;
    const auto &f_conductivity_2 = heat_face_props.medium_heat_conductivity.face_vals_axes2;
    for (auto col{0ll}; col < f_conductivity_1.cols(); ++col)
    {
      for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_1(row, col),
                   WithinRel(
                       1.0 /
                           ((grid_z.dual_nodes(row + 1ll) - grid_z.mesh_nodes(row)) /
                                heat_props.medium_heat_conductivity_axes1.value(row, col) +
                            (grid_z.mesh_nodes(row + 1ll) - grid_z.dual_nodes(row + 1ll)) /
                                heat_props.medium_heat_conductivity_axes1.value(row + 1ll, col)),
                       tol));
      }
    }

    { // col == 0
      const ptrdiff_t col = 0ll;
      for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
      {
        INFO("mesh_node: " << grid_r.mesh_nodes(col + 1ll) << ", dual_node: " << grid_r.dual_nodes(col + 1ll) << ", cement_conductivity: " << Sandface.heat_conductivity);
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       Sandface.heat_conductivity /
                           std::log(Sandface.outer_radius / Sandface.inner_radius) *
                           std::log(grid_r.dual_nodes(2ll) / grid_r.mesh_nodes(1ll)) /
                           std::log(grid_r.mesh_nodes(1ll) / grid_r.dual_nodes(1ll)),
                       tol));
      }
    }
    { // col == 1
      const ptrdiff_t col = 1ll;
      for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       1.0 /
                           (std::log(Sandface.outer_radius / Sandface.inner_radius) / Sandface.heat_conductivity +
                            std::log(grid_r.mesh_nodes(col + 1ll) / grid_r.dual_nodes(col + 1ll)) / heat_props.medium_heat_conductivity_axes1.value(row, 2ll)),
                       tol));
      }
    }
    for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
    {
      for (auto col{2ll}; col < f_conductivity_2.cols(); ++col)
      {
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       heat_props.medium_heat_conductivity_axes2.value(row, col) /
                           std::log(grid_r.dual_nodes(col + 1ll) / grid_r.dual_nodes(col)),
                       tol));
        CHECK_THAT(heat_props.medium_heat_conductivity_axes2.value(row, col),
                   WithinRel(
                       heat_props.medium_heat_conductivity_axes2.value(row, col + 1ll),
                       tol));
      }
    }
  }

  heat_face_props.apply_well(completion, well);

  // CHECK heat_face_props --- after "apply_well"
  {
    const RealType tol = 1e-12;
    // volumetric_heat_capacity
    //  const auto &capacity = heat_props.medium_vol_heatcapacity.values();
    //  const auto &porosity = core_data.porosity.log_vals;
    const auto &heat_conductivity = heat_logs.medium_heat_conductivity.log_vals;
    const auto &f_conductivity_1 = heat_face_props.medium_heat_conductivity.face_vals_axes1;
    const auto &f_conductivity_2 = heat_face_props.medium_heat_conductivity.face_vals_axes2;

    {
      const auto col{0ll}; // flow in the tube
      for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       1.0 /
                           (std::log(Tube.outer_radius / Tube.inner_radius) / Tube.heat_conductivity +
                            std::log(Annulus.outer_radius / Annulus.inner_radius) / Annulus.heat_conductivity +
                            std::log(Column.outer_radius / Column.inner_radius) / Column.heat_conductivity),
                       tol));
      }
    }

    {
      const auto col{1ll}; // flow in the tube
      for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       1.0 /
                           (std::log(Sandface.outer_radius / Sandface.inner_radius) / Sandface.heat_conductivity +
                            std::log(grid_r.mesh_nodes(col + 1ll) / grid_r.dual_nodes(col + 1ll)) / heat_conductivity(row)),
                       tol));
      }
    }

    for (auto col{2ll}; col < f_conductivity_2.cols(); ++col)
    {
      for (auto row{0ll}; row < heat_conductivity.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       heat_conductivity(row) /
                           std::log(grid_r.dual_nodes(col + 1ll) /
                                    grid_r.dual_nodes(col)),
                       tol));
      }
    }

    {
      const auto col{0ll}; // flow in the tube
      for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_1(row, col),
                   WithinRel(
                       water.heat_conductivity / (grid_z.mesh_nodes(row + 1ll) - grid_z.mesh_nodes(row)),
                       tol));
      }
    }

    {
      const auto col{1ll}; // flow in the tube
      for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_1(row, col),
                   WithinRel(
                       casing_vert_cond / (grid_z.mesh_nodes(row + 1ll) - grid_z.mesh_nodes(row)),
                       tol));
      }
    }

    for (auto col{2ll}; col < f_conductivity_1.cols(); ++col)
    {
      for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_1(row, col),
                   WithinRel(
                       1.0 / ((grid_z.dual_nodes(row + 1ll) - grid_z.mesh_nodes(row)) /
                                  heat_logs.medium_heat_conductivity(row) +
                              (grid_z.mesh_nodes(row + 1ll) - grid_z.dual_nodes(row + 1ll)) /
                                  heat_logs.medium_heat_conductivity(row + 1ll)),
                       tol));
      }
    }
  }
}
