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
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/History/RatesFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Completion.hpp>
#include <Injector/Model/ExtrudedCasingFactory.hpp>

#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/FlowField.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/make_r_stencils.hpp"
#include "includes/make_history.hpp"
#include "includes/make_geotherma.hpp"
#include "includes/get_completion.hpp"
#include "includes/generate_stencils_and_steps.hpp"
#include "includes/set_is_permeable_stencils.hpp"

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
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::FullImplicit;

/// @brief Initial temperature is assumed to be constant
struct FunctorIC : public InitialConditions::ICFunctorBase
{
  FunctorIC(const Logs::Geotherma &geotherma)
      : geotherma{geotherma}
  {
  }
  RealType operator()(const ptrdiff_t z_id, const ptrdiff_t, RealType) const override
  {
    return geotherma(z_id);
  }

protected:
  const Logs::Geotherma &geotherma;
};

template <typename Grid_t_ptr>
auto ICFactory(RealType t0, const Grid_t_ptr grid, const Logs::Geotherma &geotherma)
{
  return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{geotherma}, t0)};
}

template <typename Well_t>
struct FunctorBC : public BoundaryConditions::BCFunctorBase
{
  using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
  using ConvectionFieldFactory_t =
      GPN::FaceProperties::RatesFactory<
          Grid2D_t, Well_t, PhasePropertiesJT>;
  FunctorBC(
      const Logs::IsPermeable &is_permeable,
      const ConvectionFieldFactory_t &flow_field, // volumetric heat flow rate
      const cptr<const Grid2D_t> grid_ptr)
      : flow_field{flow_field},
        is_permeable{is_permeable},
        grid_ptr{grid_ptr}
  {
  }

  RealType operator()(const ptrdiff_t z_id, const RealType r, const RealType t) const override
  {
    if (r == grid_ptr->second_coord.dual_front())
      return flow_field.get_flow_in_axes2()(z_id, 0ll) * flow_field.get_temperature();

    if (r == grid_ptr->second_coord.dual_back())
      return 0.0;

    assert(false);
    return 0.0;
  }

  RealType operator()(const RealType z, const ptrdiff_t r_id, const RealType t) const override
  {
    if (z == grid_ptr->first_coord.dual_front())
      return flow_field.get_flow_in_axes1()(0ll, r_id) * flow_field.get_temperature();

    if (z == grid_ptr->first_coord.dual_back())
      return 0.0;

    assert(false);
    return 0.0;
  }

protected:
  const Logs::IsPermeable &is_permeable;
  const cptr<const Grid2D_t> grid_ptr;
  const ConvectionFieldFactory_t &flow_field;
};

TEST_CASE("Solver", "SelfSimilarCyl")
{
  ifstream f("heatflow_test_data.json");
  REQUIRE(f.is_open());
  json data = json::parse(f);

  /*START*/
  // input parameters
  /*fluid*/
  RealType
      viscosity{data["fluid"]["viscosity"]},
      density{data["fluid"]["density"]},
      capacity{data["fluid"]["specific_heat_capacity"]},
      heat_conductivity{data["fluid"]["heat_conductivity"]},
      joule_thomson{data["fluid"]["joule_thomson"]};
  /*collector*/
  const VR thickness{data["collector"]["thickness"].get<VR>()};
  // const ptrdiff_t nLayers{thickness.size()};
  // hydrodynamic logs
  const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
  const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"].get<VR>())};
  const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>(), 1e-12)};
  const auto RFP_weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};
  const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
  const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};
  const auto is_permeable_stencils{set_is_permeable_stencils(is_perforated_stencils, to_layers)};
  // heat logs
  const auto solid_heatconductivity_stencils{transfer_to_eigen(data["collector"]["heatConductivity"].get<VR>())};
  const auto solid_density_stencils{transfer_to_eigen(data["collector"]["solidDensity"].get<VR>())};
  const auto solid_specific_heatcapacity_stencils{transfer_to_eigen(data["collector"]["solidSpecificHeatCapacity"].get<VR>())};
  /*grid*/
  const auto
      z_minor_step{data["grid"]["z_minor_step"].get<RealType>()}; // m
  //  const ptrdiff_t rNodes{data["grid"]["rNodes"]};
  /*history*/
  const auto t_minor_step{data["history"]["t_minor_step"].get<RealType>()};
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
  const auto &grid_z{grid2D->first_coord};
  const auto &grid_r{grid2D->second_coord};

  const auto rMin{grid_r.dual_front()};
  const auto rMax{grid_r.dual_back()};

  const ExtrudedCasing extr_completion{
      VarExtrudedCasingFactory::create(completion)};

  cout << "radial dual grid stencils:\n"
       << grid_r.dual_nodes.transpose() << endl;

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

  const Logs::Rocks::CoreSampleLogs core_data{
      is_permeable_stencils,
      is_perforated_stencils,
      porosity_stencils,
      permeability_stencils,
      grid_z};

  // make fluid
  const PhasePropertiesJT water{
      FluidFactory::create_water_JT(
          Viscosity{viscosity},
          GPN::Density{density},
          GPN::SpecificHeatCapacity{capacity},
          GPN::HeatConductivity{heat_conductivity},
          JouleThomson{joule_thomson})};
  // well
  // const Well_KH well{
  //     water, core_data.is_permeable, core_data.is_perforated, core_data.permeability, well_holes, rMax};
  const auto RFP_weights{
      RFPFactory::create_from_container(
          RFP_weights_stencils,
          core_data.is_permeable)};
  const CrossFlows cross_flows{
      RFP_weights, from_coords, to_layers};
  const auto WFP_weights{
      create_WFP(
          core_data.is_perforated,
          RFP_weights,
          cross_flows)};

  const Well_Explicit well_explicit{
      core_data.is_permeable, core_data.is_perforated, RFP_weights};

  const Well_CrossFlow well{RFP_weights, WFP_weights, cross_flows};

  const Logs::Rocks::HeatLogs heat_logs{
      solid_density_stencils,
      solid_specific_heatcapacity_stencils,
      solid_heatconductivity_stencils,
      porosity_stencils,
      water,
      grid2D->first_coord};

  std::unique_ptr<const Logs::Geotherma> geotherma{
      make_unique<Logs::Geotherma>(
          make_geotherma(data, grid2D))};

  Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  // properties of material that fills the well up to the sandface
  heat_props.apply_well(extr_completion, well);

  FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};
  heat_face_props.apply_well(extr_completion, well);
  // history
  const History history{make_history(data)};
  // rates field factory
  FaceProperties::RatesFactory rates_factory{
      grid2D, well, history, water};
  // initial condition
  const auto initial_state{ICFactory(start_time, grid2D, *geotherma)};
  // boundary conditions
  const GPN::BoundaryConditions::BoundaryConditions bc{
      *grid2D,
      std::make_shared<FunctorBC<std::remove_const<decltype(well)>::type>>(
          core_data.is_permeable, rates_factory, grid2D),
      BoundaryConditions::BoundaryCondition::second};
  // solver

  using Solver_t = decltype(Solver{
      heat_face_props.medium_heat_conductivity,
      grid2D,
      heat_props.medium_vol_heatcapacity,
      rates_factory, initial_state,
      bc, start_time});

  auto solver_ptr = std::make_shared<Solver_t>(
      heat_face_props.medium_heat_conductivity,
      grid2D,
      heat_props.medium_vol_heatcapacity,
      rates_factory, initial_state,
      bc, start_time);

  const auto &solver{*solver_ptr};

  SolverManager solver_manager{history, solver_ptr};

  solver_manager.run(t_minor_step);

  // assert solution
  const double tol = 1E-11;
  const auto precision{1e-5};

  {
    string path{std::string{"flow_field.txt"}};
    ofstream f{path};
    f << (rates_factory.get_flow_in_axes1() / precision).round() * precision << endl
      << endl;
    f << (rates_factory.get_flow_in_axes2() / precision).round() * precision << endl
      << endl;
    f.close();
  }
#pragma region CHECKS
  const auto rate{rates_factory.get_rate()};
  const auto pressure{rates_factory.get_pressure()};
  // verify flow field
  const auto &v1 = rates_factory.get_flow_in_axes1(); // verticle flow
  for (auto row{0ll}; row < v1.rows(); ++row)
  {
    CHECK(v1(row, 0ll) * rate >= 0.0); // flow in tube
    CHECK(v1(row, 1ll) == 0.0);        // flow in tube-wall::annulus::column-wall
    // CHECK(col == 2ll) // flow in cement is pretty much complex
    for (auto col{3ll}; col < v1.cols(); ++col)
      CHECK(v1(row, col) == 0.0);
  }
  const auto &v2 = rates_factory.get_flow_in_axes2(); // horizontal flow
  // boundary conditions at r = 0.0
  for (auto row{0ll}, col{0ll}; row < v2.rows(); ++row)
  {
    CHECK(v2(row, col) == 0.0);
    CHECK_THAT(v1(row, col), WithinRel(v1(row + 1, col) + v2(row, col + 1ll), tol));
  }
  // compare against WFP
  const auto WFP{(water.volumetric_heat_capacity*well.get_WFP(rates_factory.get_history_record())).eval()};
  for (auto row{0ll}, col{1ll}; row < v2.rows(); ++row)
  {
    CHECK(v2(row, col) == v2(row, 2ll));
    CHECK(v2(row, col) == WFP(row));
  }
  // compare against RFP
  const auto RFP{(water.volumetric_heat_capacity*well.get_RFP(rates_factory.get_history_record())).eval()};
  for (auto row{0ll}; row < v2.rows(); ++row)
  {
    CHECK(v2(row, 3ll) == RFP(row));
    for (auto col{4ll}; col < v2.cols(); ++col)
      CHECK(v2(row, col) == v2(row, 3ll));
  }

  // flow volume balance
  for (auto row{0ll}; row < grid2D->first_coord.mesh_size(); ++row)
  {
    for (auto col{0ll}; col < grid2D->second_coord.mesh_size(); ++col)
    {
      INFO("" << "col: " << col << ", row: " << row << ", bottom: " << -v1(row + 1ll, col) << ", top: " << v1(row, col) << ", right: " << v2(row, col + 1ll) << ", left: " << -v2(row, col));
      CHECK_THAT(-v1(row + 1ll, col) + v1(row, col), WithinRel(v2(row, col + 1ll) - v2(row, col), tol));
    }
  }
  // maximum principle
  const auto &[times, states] = solver.solution();
  // for (size_t i{0ll}; i < times.size(); ++i)
  // {
  //   const auto &state = states[i];
  //   for (auto row{0ll}; row < state.rows(); ++row)
  //   {
  //     CHECK(state(row, 0ll) >= inlet_temperature - tol);
  //     for (auto col{1ll}; col < state.cols(); ++col)
  //     {
  //       INFO("time: " << i << ", col: " << col << ", row: " << row);
  //       CHECK(state(row, col) >= inlet_temperature);
  //     }
  //   }
  // }
  // for (size_t i{1ull}; i < times.size(); ++i)
  // {
  //   const auto &state = states[i];
  //   for (auto row{0ll}; row < state.rows(); ++row)
  //   {
  //     for (auto col{1ll}; col < state.cols(); ++col)
  //     {
  //       CHECK((states[i](row, col) - states[i - 1ull](row, col)) / (states[i](row, col) + states[i - 1ull](row, col)) <= tol);
  //     }
  //   }
  // }

  // overall heat balance
  RealType cur_heat_incr = 0.0;
  RealType cum_inlet_heat = 0.0;
  //  cout << "volumetric heat capacity\n"
  //       << heat_props.medium_vol_heatcapacity.its_values << endl;

  for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
  {
    cur_heat_incr +=
        ((states[t].cur_state - states[t - 1ll].cur_state) *
         heat_props.medium_vol_heatcapacity.its_values * grid2D->volumes())
            .sum();
    cum_inlet_heat +=
        (times[t] - times[t - 1ll]) *
        history.rates(t - 1ll) *
        water.volumetric_heat_capacity * (history.temps(t - 1ll) /*- initial_temperature*/);

    RealType rel_tol = std::abs(2.0 * (cur_heat_incr - cum_inlet_heat) / (cur_heat_incr + cum_inlet_heat));
    //    CHECK(rel_tol < 0.05);
  }
#pragma endregion
  {
    const auto &state = states.back();
    std::string path{std::string{"T_"} + std::to_string(0) + std::string{".txt"}};
    std::ofstream f{path};

    f << ((state.cur_state /*- initial_temperature*/) / precision).round() * precision;
    f.close();
  }
}