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
#include <Injector/Model/Completion.hpp>
#include <Injector/Model/ExtrudedCasingFactory.hpp>

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

TEST_CASE("apply_well_test", "apply_well_test")
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
  REQUIRE(heat_conductivity > 0.0);
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
  for (const auto &r : completion.sandwich)
    REQUIRE(r.heat_conductivity > 0.0);
  /*END*/

  // make grid2D
  // r_stencils
  const WellHoles well_holes{WellHolesFactory::create(completion)};
  const VR r_stencils{make_r_stencils(data, well_holes)};
  const RealType &rMax = r_stencils.back();
  const RealType &rMin = r_stencils.front();

  for (const auto &r : completion.sandwich)
  {
    CHECK_THAT(r.linear_heat_capacity,
               WithinRel(
                   r.area() * r.volumetric_heat_capacity,
                   tol));
  }

  // CHECK r_stencils
  {
    CHECK(r_stencils[0ull] == 0.0);
    CHECK(r_stencils[1ull] == completion.flow_radius);
    CHECK(r_stencils[2ull] == completion.column_outer_radius);
    CHECK(r_stencils[3ull] == completion.sandface_radius);
    CHECK(r_stencils[4ull] > r_stencils[3ll]);
  }

  // z-refiner
  RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
  // the grid itself
  const auto grid2D{
      Grids::CylinderGridFactory::create(
          refiner,
          Grids::Factory::generate_dual_grid_stencils_from_steps(
              0.0, thickness),
          r_stencils)};
  const auto &grid_z{grid2D->first_coord};
  const auto &grid_r{grid2D->second_coord};

  const ExtrudedCasing extr_completion{ExtrudedCasingFactory::create(completion, grid_z)};
  const auto &Flow{extr_completion[MaterialType::Flow]};
  const auto &Tube{extr_completion[MaterialType::Tube]};
  const auto &Annulus{extr_completion[MaterialType::Annulus]};
  const auto &Column{extr_completion[MaterialType::Column]};
  const auto &Sandface{extr_completion[MaterialType::Cement]};
  for (auto i{0ll}; i < grid_z.mesh_size(); ++i)
  {
    CHECK_THAT(grid2D->face_area_axes2(i),
               WithinRel(2 * numbers::pi *
                             (grid_z.dual_nodes(i + 1ll) - grid_z.dual_nodes(i)),
                         tol));
    CHECK(Flow.thickness(i) > 0.0);
    CHECK(Column.thickness(i) > 0.0);
    CHECK(Sandface.thickness(i) > 0.0);
    CHECK(Flow.outer_radius(i) == Tube.inner_radius(i));
    CHECK(Tube.outer_radius(i) == Annulus.inner_radius(i));
    CHECK(Annulus.outer_radius(i) == Column.inner_radius(i));
    CHECK(Column.outer_radius(i) == Sandface.inner_radius(i));
    if (grid_z.mesh_nodes(i) < completion[MaterialType::Tube].depth)
    {
      CHECK(Tube.thickness(i) > 0.0);
      CHECK(Annulus.thickness(i) > 0.0);
    }
    else
    {
      CHECK(Tube.thickness(i) == 0.0);
      CHECK(Annulus.thickness(i) == 0.0);
    }
  }

  const Eigen::ArrayX<RealType> casing_vert_cond{
      std::numbers::pi *
      (Tube.heat_conductivity * Tube.thickness * (Tube.inner_radius + Tube.outer_radius) +
       Annulus.heat_conductivity * Annulus.thickness * (Annulus.inner_radius + Annulus.outer_radius) +
       Column.heat_conductivity * Column.thickness * (Column.inner_radius + Column.outer_radius)) /
      (std::numbers::pi *
       (Column.outer_radius + Tube.inner_radius) *
       (Column.outer_radius - Tube.inner_radius))};
  const Eigen::ArrayX<RealType> cement_vert_cond{
      std::numbers::pi * Sandface.thickness * (Sandface.inner_radius + Sandface.outer_radius) * Sandface.heat_conductivity /
      (std::numbers::pi * Sandface.thickness * (Sandface.inner_radius + Sandface.outer_radius))};
  {
    const Eigen::ArrayX<RealType> cement_temp{extr_completion.integral_vertical_cement_heat_conductivity()};
    const Eigen::ArrayX<RealType> casing_temp{extr_completion.integral_vertical_casing_heat_conductivity()};
    const Eigen::ArrayX<RealType> casing_area{extr_completion.casing_area()};
    const Eigen::ArrayX<RealType> cement_area{extr_completion.cement_area()};

    for (auto i{0ll}; i < casing_vert_cond.size(); ++i)
    {
      // check verticle conductivity
      CHECK_THAT(casing_vert_cond(i), WithinRel(casing_temp(i), tol));
      CHECK_THAT(cement_vert_cond(i), WithinRel(cement_temp(i), tol));
      // check const verticle conductivity in cement
      CHECK_THAT(cement_vert_cond(0ll), WithinRel(cement_vert_cond(i), tol));
      // check area
      CHECK_THAT((std::numbers::pi *
                  (Column.outer_radius(i) + Tube.inner_radius(i)) * (Column.outer_radius(i) - Tube.inner_radius(i))),
                 WithinRel(casing_area(i), tol));
      // check area
      CHECK_THAT((std::numbers::pi *
                  (Sandface.outer_radius(i) + Sandface.inner_radius(i)) * Sandface.thickness(i)),
                 WithinRel(cement_area(i), tol));
    }
  }

  const Eigen::ArrayX<RealType> cement_heat_cap{
      completion[MaterialType::Cement].volumetric_heat_capacity * Sandface.area()};
  {
    const auto &cement_temp{Sandface.linear_heat_capacity};
    for (auto i{0ll}; i < casing_vert_cond.size(); ++i)
    {
      // check linear heat capacity
      CHECK_THAT(cement_heat_cap(i), WithinRel(cement_temp(i), tol));
      // check const linear heat capacity in cement
      CHECK_THAT(cement_heat_cap(0ll), WithinRel(cement_heat_cap(i), tol));
    }
  }

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
          GPN::Density{density},
          GPN::SpecificHeatCapacity{capacity},
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
      // weight value
      CHECK_THAT(w(row), WithinRel(weights_stencils(idx), tol));
    }
    CHECK_THAT(w.sum(), WithinRel(1.0, tol));
    CHECK_THAT(weights_stencils.sum(), WithinRel(1.0, tol));
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
  heat_props.apply_well(extr_completion, well);

  // CHECK heat_props --- after "apply_well"
  {
    const auto &capacity = heat_props.medium_vol_heatcapacity.values();
    const auto &heat_conductivity_1 = heat_props.medium_heat_conductivity_axes1.values();
    const auto &heat_conductivity_2 = heat_props.medium_heat_conductivity_axes2.values();
    const auto flow_area{Flow.area()};

    const Eigen::ArrayX<RealType> temp{Flow.volumetric_heat_capacity *
                                       flow_area *
                                       (grid_z.dual_nodes.tail(capacity.rows()) - grid_z.dual_nodes.head(capacity.rows()))};
    const Eigen::ArrayX<RealType> flow_face_ratio{extr_completion.flow().area() /
                                                  grid2D->face_area_axes1(0ll)};
    const Eigen::ArrayX<RealType> casing_face_ratio{extr_completion.casing_area() /
                                                    grid2D->face_area_axes1(1ll)};
    const Eigen::ArrayX<RealType> cement_face_ratio{extr_completion.cement_area() /
                                                    grid2D->face_area_axes1(2ll)};

    for (auto row{0ll}; row < capacity.rows(); ++row)
    {
      { // col == 0 --- fluid flow
        const auto col{0ll};
        // flow heat capacity is equal to fluid-water heat capacity
        CHECK_THAT(capacity(row, col) * grid2D->volume(row, col),
                   WithinRel(
                       Flow.linear_heat_capacity(row) *
                           (grid_z.dual_nodes(row + 1ll) - grid_z.dual_nodes(row)),
                       tol));
        CHECK_THAT(capacity(row, col),
                   WithinRel(
                       extr_completion.flow().volumetric_heat_capacity * flow_face_ratio(row), tol));
        // vertical heat conductivity is equal to water corrected for the face ratio
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       water.heat_conductivity *
                           flow_area(row) / grid2D->face_area_axes1(col),
                       tol));
        // radial heat conductivity of flowing water is infinity
        CHECK(std::isinf(heat_conductivity_2(row, col)));
      }
      { // col == 1 --- casing
        const auto col{1ll};
        CHECK_THAT(capacity(row, col) * grid2D->volume(row, col),
                   WithinRel(
                       (grid_z.dual_nodes(row + 1ll) - grid_z.dual_nodes(row)) *
                           (Tube.linear_heat_capacity(row) +
                            Annulus.linear_heat_capacity(row) +
                            Column.linear_heat_capacity(row)),
                       tol));
        // vertical heat conductivity in the casing corrected for the face ratio
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       casing_vert_cond(row) * casing_face_ratio(row),
                       tol));
        // radial heat conductivity in the casing is undefined
      }
      { // col == 2 --- cement
        const auto col{2ll};
        CHECK_THAT(capacity(row, col) * grid2D->volume(row, col),
                   WithinRel(
                       (grid_z.dual_nodes(row + 1ll) - grid_z.dual_nodes(row)) *
                           Sandface.linear_heat_capacity(row),
                       tol));
        CHECK_THAT(heat_conductivity_1(row, col),
                   WithinRel(
                       cement_vert_cond(row) * cement_face_ratio(row),
                       tol));
        CHECK_THAT(heat_conductivity_2(row, col),
                   WithinRel(
                       Sandface.heat_conductivity,
                       tol));
      }

      for (auto col{3ll}; col < capacity.cols(); ++col)
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
  heat_face_props.apply_well(extr_completion, well);

  // CHECK heat_face_props --- after "apply_well"
  {
    const RealType tol = 1e-12;
    // volumetric_heat_capacity
    //  const auto &capacity = heat_props.medium_vol_heatcapacity.values();
    //  const auto &porosity = core_data.porosity.log_vals;
    const auto &heat_conductivity = heat_logs.medium_heat_conductivity.log_vals;
    const auto &f_conductivity_1 = heat_face_props.medium_heat_conductivity.face_vals_axes1;
    const auto &f_conductivity_2 = heat_face_props.medium_heat_conductivity.face_vals_axes2;

    // CHECK f_conductivity_2
    {
      // corrected position of r = r_1 -- inside the casing-sandwich
      const auto r1{extr_completion.radial_node_position()};
      CHECK(r1.rows() == Tube.inner_radius.rows());
      CHECK(r1.rows() == Annulus.inner_radius.rows());
      CHECK(r1.rows() == Column.inner_radius.rows());
      CHECK(r1.rows() == Sandface.inner_radius.rows());
      // heat resistivity at the face between the flow and the casing-sandwich
      Eigen::ArrayX<RealType> zeta_0(r1.size());
      for (auto id{0ll}; id < zeta_0.size(); ++id)
      {
        CHECK(r1(id) > Tube.inner_radius(id));
        CHECK(r1(id) < Column.outer_radius(id));
        if (completion[MaterialType::Tube].depth < grid_z.mesh_nodes(id))
        {
          CHECK(r1(id) > Column.inner_radius(id));
          CHECK(r1(id) < Column.outer_radius(id));
        }
        zeta_0(id) =
            r1(id) < Tube.outer_radius(id)
                ? std::log(r1(id) / Tube.inner_radius(id)) / Tube.heat_conductivity
            : r1(id) < Annulus.outer_radius(id)
                ? 1 / Tube.radial_heat_conductivity(id) + std::log(r1(id) / Tube.outer_radius(id)) / Annulus.heat_conductivity
                // r1(id) < Column.outer-radius(id)
                : 1 / Tube.radial_heat_conductivity(id) + 1 / Annulus.radial_heat_conductivity(id) + std::log(r1(id) / Column.inner_radius(id)) / Column.heat_conductivity;
      }

      {
        const auto col{0ll}; // face between flow and tube wall
        for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
        {
          CHECK_THAT(f_conductivity_2(row, col),
                     WithinRel(1.0 / zeta_0(row),
                               tol));
        }
      }

      // heat resistivity at the face between the casing-sandwich and the cement
      const Eigen::ArrayX<RealType> zeta_1{
          1 / extr_completion.integral_casing_radial_heat_conductivity() +
          log(grid_r.mesh_nodes(2ll) / Column.outer_radius) / Sandface.heat_conductivity -
          zeta_0};
      {
        const auto col{1ll}; // flow in the tube
        for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
        {
          CHECK_THAT(f_conductivity_2(row, col),
                     WithinRel(
                         1.0 / zeta_1(row),
                         tol));
        }
      }

      {
        const auto col{2ll}; // face between cement and rocks
        for (auto row{0ll}; row < f_conductivity_2.rows(); ++row)
        {
          const auto zeta_2{
              log(Sandface.outer_radius(row) / grid_r.mesh_nodes(col)) / Sandface.heat_conductivity +
              log(grid_r.mesh_nodes(col + 1ll) / Sandface.outer_radius(row)) / heat_conductivity(row)};
          CHECK_THAT(f_conductivity_2(row, col),
                     WithinRel(
                         1.0 / zeta_2,
                         tol));
        }
      }

      for (auto col{3ll}; col < f_conductivity_2.cols(); ++col)
      {
        for (auto row{0ll}; row < heat_conductivity.rows(); ++row)
        {
          CHECK_THAT(f_conductivity_2(row, col),
                     WithinRel(
                         heat_conductivity(row) /
                             std::log(grid_r.mesh_nodes(col + 1ll) /
                                      grid_r.mesh_nodes(col)),
                         tol));
        }
      }
    }
    // CHECK f_conductivity_1
    {
      {
        const auto col{0ll}; // flow in the tube
        const Eigen::ArrayX<RealType> flow_area{Flow.area()};
        for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
        {
          CHECK_THAT(f_conductivity_1(row, col) * grid2D->face_area_axes1(col),
                     WithinRel(
                         1.0 / ((grid_z.dual_nodes(row + 1ll) - grid_z.mesh_nodes(row)) /
                                    (Flow.heat_conductivity * flow_area(row)) +
                                (grid_z.mesh_nodes(row + 1ll) - grid_z.dual_nodes(row + 1ll)) /
                                    (Flow.heat_conductivity * flow_area(row + 1ll))),
                         tol));
        }
      }

      {
        const auto col{1ll}; // conductivity along the casing
        const Eigen::ArrayX<RealType> temp{
            extr_completion.integral_vertical_casing_heat_conductivity() *
            extr_completion.casing_area()};
        for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
        {
          CHECK_THAT(f_conductivity_1(row, col) * grid2D->face_area_axes1(col),
                     WithinRel(
                         1.0 /
                             ((grid_z.dual_nodes(row + 1ll) - grid_z.mesh_nodes(row)) /
                                  temp(row) +
                              (grid_z.mesh_nodes(row + 1ll) - grid_z.dual_nodes(row + 1ll)) /
                                  temp(row + 1ll)),
                         tol));
        }
      }

      {
        const auto col{2ll}; // conductivity along the cement
        const Eigen::ArrayX<RealType> temp{
            extr_completion.integral_vertical_cement_heat_conductivity() *
            extr_completion.cement_area()};
        for (auto row{0ll}; row < f_conductivity_1.rows(); ++row)
        {
          CHECK_THAT(f_conductivity_1(row, col),
                     WithinRel(
                         cement_vert_cond(row) / (grid_z.mesh_nodes(row + 1ll) - grid_z.mesh_nodes(row)),
                         tol));
          CHECK_THAT(f_conductivity_1(row, col) * grid2D->face_area_axes1(col),
                     WithinRel(
                         1.0 /
                             ((grid_z.dual_nodes(row + 1ll) - grid_z.mesh_nodes(row)) /
                                  temp(row) +
                              (grid_z.mesh_nodes(row + 1ll) - grid_z.dual_nodes(row + 1ll)) /
                                  temp(row + 1ll)),
                         tol));
        }
      }

      for (auto col{3ll}; col < f_conductivity_1.cols(); ++col)
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
}
