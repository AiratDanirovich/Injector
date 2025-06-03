#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <numbers>
#include <cmath>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/History/History.hpp>
#include <Injector/History/RatesFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include <Injector/Properties/FieldsFactory.hpp>

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;
using json = nlohmann::json;

using namespace std;
using namespace GPN;
using namespace GPN::Logs;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

/// @brief Initial temperature is assumed to be constant
struct FunctorIC
{
  FunctorIC(const RealType val)
      : val{val}
  {
  }
  RealType operator()(RealType z, RealType r, RealType t0) const
  {
    return val;
  }

protected:
  const RealType val;
};

template <typename Grid_t_ptr>
auto ICFactory(RealType t0, const Grid_t_ptr grid, const RealType val)
{
  return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{val}, t0)};
}

struct FunctorBC : public BoundaryConditions::BCFunctorBase
{
  using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
  using ConvectionFieldFactory_t =
      GPN::FaceProperties::RatesFactory<
          Grid2D_t, Well_KH, PhaseProperties>;
  FunctorBC(
      const RealType inlet_temp,
      const Logs::IsPermeable &is_permeable,
      const ConvectionFieldFactory_t &flow_field, // volumetric heat flow rate
      const cptr<const Grid2D_t> grid_ptr)
      : inlet_temp{inlet_temp},
        flow_field{flow_field},
        is_permeable{is_permeable},
        grid_ptr{grid_ptr}
  {
  }

  RealType operator()(ptrdiff_t z_id, RealType r, RealType t) const override
  {
    if (r == grid_ptr->second_coord.dual_front())
      return flow_field.get_flow_in_axes2()(z_id, 0ll) * inlet_temp;

    if (r == grid_ptr->second_coord.dual_back())
      return 0.0;

    assert(false);
    return 0.0;
  }

  RealType operator()(RealType z, ptrdiff_t r_id, RealType t) const override
  {
    if (z == grid_ptr->first_coord.dual_front())
      return flow_field.get_flow_in_axes1()(0ll, r_id) * inlet_temp;

    if (z == grid_ptr->first_coord.dual_back())
      return 0.0;

    assert(false);
    return 0.0;
  }

protected:
  RealType inlet_temp;
  const Logs::IsPermeable &is_permeable;
  const cptr<const Grid2D_t> grid_ptr;
  const ConvectionFieldFactory_t &flow_field;
};

using VR = std::vector<RealType>;
LogValuesContainer transfer_to_eigen(const VR &data, const RealType factor = 1.0)
{
  LogValuesContainer out(data.size());
  for (auto i{0ull}; i < data.size(); ++i)
    out(i) = factor * data[i];
  return out;
}

VR generate_stencils(RealType t0, RealType t1, RealType t_step_major)
{
  auto segm_count{static_cast<size_t>(std::ceil(t1 - t0) / t_step_major)};
  double step = (t1 - t0) / (segm_count);
  VR out(segm_count + 1ll);

  for (auto i{0ull}; i < out.size(); ++i)
    out[i] = t0 + i * step;
  return out;
}

VR generate_steps(const VR &dual_nodes)
{
  VR out(dual_nodes.size() - 1ll);

  for (auto i{0ull}; i < out.size(); ++i)
    out[i] = dual_nodes[i + 1] - dual_nodes[i];
  return out;
}

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
      capacity{data["fluid"]["specificHeatCapacity"]},
      heat_conductivity{data["fluid"]["heatConductivity"]};
  /*collector*/
  const VR thickness = data["collector"]["thickness"];
  //  const ptrdiff_t nLayers{thickness.size()};
  // hydrodynamic logs
  const auto is_permeable_stencils{transfer_to_eigen(data["collector"]["is_permeable"])};
  const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"])};
  const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"], 1e-12)};
  // heat logs
  const VR heatconductivity_stencils = data["collector"]["heatConductivity"];
  const auto solid_density_stencils{transfer_to_eigen(data["collector"]["solidDensity"])};
  const auto solid_specific_heatcapacity_stencils{transfer_to_eigen(data["collector"]["solidSpecificHeatCapacity"])};
  /*grid*/
  const RealType
      rMin{data["grid"]["r_start"]},
      rMax{data["grid"]["r_end"]},
      zTop{data["grid"]["ztop"]}; // m
  const ptrdiff_t rNodes{data["grid"]["rNodes"]};
  /*history*/
  const RealType
      t0{data["history"]["t_start"]},
      t1{data["history"]["t_end"]};
  RealType t_major_step{data["history"]["t_major_step"]};
  RealType t_minor_step{data["history"]["t_minor_step"]};
  t_major_step = std::min(t1 - t0, t_major_step);
  t_minor_step = std::min(t_minor_step, t_major_step);
  const VR t_stencils{generate_stencils(t0, t1, t_major_step)};
  /*temperatures*/
  const RealType well_rate{data["history"]["wellRate"]}; // m^3/s
  const RealType initial_temperature{data["collector"]["initTemperature"]};
  const RealType inlet_temperature{data["history"]["inletTemperature"]};
  /*well*/
  const RealType hole_radius{data["well"]["hole_radius"]};
  /*END*/

  REQUIRE(t1 > t0);
  REQUIRE(t_minor_step <= t_major_step);
  //  REQUIRE(rMin < hole_radius);

  // make grid2D
  VR r_stencils;
  //  r_stencils.push_back(rMin);
  auto temp = Grids::Factory::generate_dual_grid_stencils_uniform(
      Segment{hole_radius, rMax}, rNodes);
  r_stencils.insert(r_stencils.end(), temp.begin(), temp.end());

  const auto grid2D{
      Grids::CylinderGridFactory::create(
          Grids::Factory::generate_dual_grid_stencils_from_steps(
              zTop, thickness),
          r_stencils)};
  const auto &grid{grid2D->first_coord};

  cout << "radial grid:\n"
       << grid2D->second_coord.dual_nodes.transpose() << endl;
  cout << "vertical grid:\n"
       << grid2D->first_coord.dual_nodes.transpose() << endl;

  cout << "radial grid cell centers:\n"
       << grid2D->second_coord.mesh_nodes.transpose() << endl;
  cout << "vertical grid cell centers:\n"
       << grid2D->first_coord.mesh_nodes.transpose() << endl;

  cout << "radial grid mesh steps:\n"
       << grid2D->second_coord.mesh_steps.transpose() << endl;
  cout << "vertical grid mesh steps:\n"
       << grid2D->first_coord.mesh_steps.transpose() << endl;

  const Logs::Rocks::CoreSampleLogs core_data{
      is_permeable_stencils,
      porosity_stencils,
      permeability_stencils,
      grid};

  // make fluid
  const PhaseProperties water{
      FluidFactory::create_water(
          Viscosity{viscosity},
          Density{density},
          SpecificHeatCapacity{capacity},
          GPN::HeatConductivity{heat_conductivity})};
  // well
  const Well_KH well{
      water, core_data.is_permeable, core_data.permeability};

  const Logs::Rocks::HeatLogs heat_logs{
      solid_density_stencils,
      solid_specific_heatcapacity_stencils,
      heatconductivity_stencils,
      porosity_stencils,
      water,
      grid2D->first_coord};
  Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  heat_props.apply_well(well, water);

  const FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};
  // history
  const std::vector<RealType> time_steps{generate_steps(t_stencils)};
  const std::vector<RealType> rates(time_steps.size(), well_rate);
  const History history{
      HistoryFactory::create(time_steps, rates)};
  // rates field factory
  FaceProperties::RatesFactory rates_factory{
      grid2D, well, history, water};
  // initial condition
  const auto initial_state{ICFactory(t0, grid2D, initial_temperature)};
  // boundary conditions
  const GPN::BoundaryConditions::BoundaryConditions bc{
      *grid2D,
      std::make_shared<FunctorBC>(
          inlet_temperature, core_data.is_permeable, rates_factory, grid2D),
      BoundaryConditions::BoundaryCondition::second};
  // solver
  using Solver_t = decltype(Solver{
      heat_face_props.heat_conductivity,
      grid2D,
      heat_props.medium_vol_heatcapacity,
      rates_factory, initial_state,
      bc, t0});

  auto solver_ptr = std::make_shared<Solver_t>(
      heat_face_props.heat_conductivity,
      grid2D,
      heat_props.medium_vol_heatcapacity,
      rates_factory, initial_state,
      bc, t0);

  const auto &solver{*solver_ptr};

  SolverManager solver_manager{history, solver_ptr};

  solver_manager.run(t_minor_step);

  // assert solution
  const double tol = 1E-12;
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

  // verify flow field
  const auto &v1 = rates_factory.get_flow_in_axes1();
  for (auto row{0ll}; row < v1.rows(); ++row)
  {
    CHECK(v1(row, 0ll) >= 0.0);
    for (auto col{1ll}; col < v1.cols(); ++col)
      CHECK(v1(row, col) == 0.0);
  }
  const auto &v2 = rates_factory.get_flow_in_axes2();
  for (auto col{2ll}; col < v2.cols(); ++col)
    for (auto row{0ll}; row < v2.rows(); ++row)
      CHECK(v2(row, col) == v2(row, 1ll));

  for (auto row{0ll}, col{0ll}; row < v2.rows(); ++row)
  {
    CHECK(v2(row, col) == 0.0);
    CHECK_THAT(v1(row, col), WithinRel(v1(row + 1, col) + v2(row, col + 1ll), tol));
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
  for (size_t i{0ll}; i < times.size(); ++i)
  {
    const auto &state = states[i];
    for (auto row{0ll}; row < state.rows(); ++row)
    {
      CHECK(state(row, 0ll) >= inlet_temperature - tol);
      for (auto col{1ll}; col < state.cols(); ++col)
      {
        INFO("time: " << i << ", col: " << col << ", row: " << row);
        CHECK(state(row, col) >= inlet_temperature);
      }
    }
  }
  for (size_t i{1ull}; i < times.size(); ++i)
  {
    const auto &state = states[i];
    for (auto row{0ll}; row < state.rows(); ++row)
    {
      for (auto col{1ll}; col < state.cols(); ++col)
      {
        CHECK(states[i](row, col) <= states[i - 1ull](row, col) + tol);
      }
    }
  }

  // overall heat balance
  RealType cur_heat_incr = 0.0;
  RealType cum_inlet_heat = 0.0;
  cout << "volumetric heat capacity\n"
       << heat_props.medium_vol_heatcapacity.its_values << endl;

  for (auto t{1ll}; t < (ptrdiff_t)times.size(); ++t)
  {
    cur_heat_incr =
        ((states[t].cur_state - states[0ll].cur_state) *
         heat_props.medium_vol_heatcapacity.its_values * grid2D->volumes())
            .sum();
    cum_inlet_heat =
        (times[t] - times[0ll]) *
        history.rates(t - 1ll) *
        water.volumetric_heat_capacity * (inlet_temperature - initial_temperature);

    RealType rel_tol = std::abs(2.0 * (cur_heat_incr - cum_inlet_heat) / (cur_heat_incr + cum_inlet_heat));
    CHECK(rel_tol < 0.1);
  }

  {
    const auto &state = states.back();
    std::string path{std::string{"T_"} + std::to_string(0) + std::string{".txt"}};
    std::ofstream f{path};

    f << ((state.cur_state - initial_temperature) / precision).round() * precision;
    f.close();
  }
}