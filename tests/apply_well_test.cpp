// #include <memory>
#include <iostream>
#include <fstream>
#include <vector>
// #include <string>
// #include <numbers>
// #include <cmath>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/Grids2D.hpp>
// #include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Model/WellFactory.hpp>
// #include <Injector/Properties/Factory.hpp>
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

  // z-refiner
  RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
  // the grid itself
  const auto grid2D{
      Grids::CylinderGridFactory::create(refiner,
                                         Grids::Factory::generate_dual_grid_stencils_from_steps(
                                             0.0, thickness),
                                         r_stencils)};
  const auto &grid{grid2D->first_coord};
  const auto &grid_r{grid2D->second_coord};

  cout << "radial dual grid stencils:\n"
       << transfer_to_eigen(r_stencils).transpose() << endl;

  const Logs::Rocks::CoreSampleLogs core_data{
      is_permeable_stencils,
      is_perforated_stencils,
      porosity_stencils,
      permeability_stencils,
      grid};

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
              grid.dual_stencils.dual_nodes,
              grid.mesh_nodes(row))};
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
          grid)};

  // CHECK weights
  {
    const auto &w = weights.log_vals;
    for (auto row{0ll}; row < w.rows(); ++row)
    {
      // find val in coarse grid
      const ptrdiff_t idx{
          get_dist(
              grid.dual_stencils.dual_nodes,
              grid.mesh_nodes(row))};
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
              grid.dual_stencils.dual_nodes,
              grid.mesh_nodes(row))};
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

  // properties of material that fills the well up to the sandface
  heat_props.apply_well(completion, well);

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
        const auto &sandface = completion[MaterialType::CementOuter];
        CHECK_THAT(capacity(row, col),
                   WithinRel(
                       (completion[MaterialType::Tube].linear_heat_capacity +
                        completion[MaterialType::Annulus].linear_heat_capacity +
                        completion[MaterialType::Column].linear_heat_capacity +
                        completion[MaterialType::CementInner].linear_heat_capacity +
                        sandface.linear_heat_capacity) /
                           completion.area(),
                       tol));
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       (completion[MaterialType::Tube].integral_vertical_heat_conductivity +
                        completion[MaterialType::Annulus].integral_vertical_heat_conductivity +
                        completion[MaterialType::Column].integral_vertical_heat_conductivity +
                        completion[MaterialType::CementInner].integral_vertical_heat_conductivity +
                        sandface.integral_vertical_heat_conductivity) /
                           completion.area(),
                       tol));
        CHECK_THAT(heat_conductivity_2(row, col),
                   WithinRel(
                       sandface.heat_conductivity /
                           std::log(sandface.outer_radius / sandface.inner_radius) *
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

  cout << endl
       << "heat_conductivity_axes1:\n"
       << heat_props.medium_heat_conductivity_axes1.values() << endl
       << endl
       << "heat_conductivity_axes2:\n"
       << heat_props.medium_heat_conductivity_axes2.values()
       << endl;

  // CHECK heat_face_props
  {
    const auto &f_conductivity_1 = heat_face_props.medium_heat_conductivity.face_vals_axes1;
    const auto &f_conductivity_2 = heat_face_props.medium_heat_conductivity.face_vals_axes2;
    for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
    {
      for (auto col{0ll}; col < f_conductivity_1.cols(); ++col)
      {
        CHECK_THAT(f_conductivity_1(row, col),
                   WithinRel(
                       1.0 /
                           ((grid.dual_nodes(row + 1ll) - grid.mesh_nodes(row)) /
                                heat_props.medium_heat_conductivity_axes1.value(row, col) +
                            (grid.mesh_nodes(row + 1ll) - grid.dual_nodes(row + 1ll)) /
                                heat_props.medium_heat_conductivity_axes1.value(row + 1ll, col)),
                       tol));
      }
    }

    { // col == 0
      const ptrdiff_t col = 0ll;
      const auto &CementOuter = completion[MaterialType::CementOuter];
      for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
      {
        INFO("mesh_node: " << grid_r.mesh_nodes(col + 1ll) << ", dual_node: " << grid_r.dual_nodes(col + 1ll) << ", cement_conductivity: " << CementOuter.heat_conductivity);
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       1.0 /
                           (std::log(grid_r.mesh_nodes(col + 1ll) / grid_r.dual_nodes(col + 1ll)) / CementOuter.heat_conductivity),
                       tol));
      }
    }
    { // col == 1
      const ptrdiff_t col = 1ll;
      const auto &CementOuter = completion[MaterialType::CementOuter];
      for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
      {
        CHECK_THAT(f_conductivity_2(row, col),
                   WithinRel(
                       1.0 /
                           (std::log(CementOuter.outer_radius / CementOuter.inner_radius) / CementOuter.heat_conductivity +
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

  // CHECK heat_props
  // {
  //   const RealType tol = 1e-12;
  //   // volumetric_heat_capacity
  //   const auto &capacity = heat_props.medium_vol_heatcapacity.values();
  //   const auto &porosity = core_data.porosity.log_vals;
  //   for (auto row{0ll}; row < capacity.rows(); ++row)
  //   {
  //     for(auto col{2ll}; col < capacity.cols(); ++ col)
  //     {
  //       CHECK_THAT(capacity(row, col),
  //       WithinRel(porosity(row)*water.volumetric_heat_capacity + (1-porosity(row))*heat_logs.solid_vol_heatcapacity(row), tol));
  //     }
  //   }
  // }
}