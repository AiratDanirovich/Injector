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
      RealType inlet_temp,
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
      return flow_field.get_flow_in_axes1()(z_id, 0ll) * inlet_temp;

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

TEST_CASE("Solver", "SelfSimilarCyl")
{
  ifstream f("../../../tests/test_data/heatflow_test_data.json");
  REQUIRE(f.is_open());
  json data = json::parse(f);

  /*START*/
  // input parameters
  /*fluid*/
  RealType
      viscosity{data["fluid"]["viscosity"]},
      density{data["fluid"]["density"]},
      capacity{data["fluid"]["specificHeatCapacity"]};
  /*collector*/
  const ptrdiff_t nLayers{data["collector"]["nLayers"]};
  const VR thickness(nLayers, data["collector"]["thickness"]);
  // hydrodynamic logs
  auto is_permeable_stencils{
      LogValuesContainer::Constant(nLayers, 1.0).eval()};
  auto porosity_stencils{
      LogValuesContainer::Constant(nLayers, data["collector"]["porosity"]).eval()};
  auto permeability_stencils{
      LogValuesContainer::Constant(nLayers, data["collector"]["permeability"]).eval()};
  // heat logs
  const VR heatconductivity_stencils(nLayers, data["collector"]["heatConductivity"]);
  const auto solid_density_stencils{
      LogValuesContainer::Constant(nLayers, data["collector"]["solidDensity"])};
  const auto solid_specific_heatcapacity_stencils{
      LogValuesContainer::Constant(nLayers, data["collector"]["solidSpecificHeatCapacity"])};
  /*grid*/
  const RealType
      rMin{data["grid"]["r_start"]},
      rMax{data["grid"]["r_end"]},
      zTop{data["grid"]["ztop"]}; // m
  const ptrdiff_t rNodes{data["grid"]["rNodes"]};
  /*history*/
  const RealType t0{data["history"]["t_start"]},
      t1{data["history"]["t_end"]};
  REQUIRE(t1 > t0);
  const ptrdiff_t time_steps_nmbr{static_cast<ptrdiff_t>(ceil(
      (t1 - t0) / (double)data["history"]["t_step"]))};
  const RealType time_step{(t1 - t0) / time_steps_nmbr};
  const VR time_intervals(time_steps_nmbr, time_step);
  const VR t_stencils(
      Grids::Factory::generate_dual_grid_stencils_from_steps(
          t0, time_intervals));
  /*temperatures*/
  const RealType well_rate{data["history"]["wellRate"]}; // m^3/s
  const RealType initial_temperature{data["collector"]["initTemperature"]};
  const RealType inlet_temperature{data["history"]["inletTemperature"]};
  /*END*/

  // make grid2D
  const auto grid2D{
      Grids::CylinderGridFactory::create(
          Grids::Factory::generate_dual_grid_stencils_from_steps(
              zTop, thickness),
          Grids::Factory::generate_dual_grid_stencils_uniform(
              Segment{rMin, rMax}, rNodes))};
  const auto &grid{grid2D->first_coord};

  is_permeable_stencils(nLayers / 4) = 0.0;
  is_permeable_stencils(nLayers / 2) = 0.0;
  is_permeable_stencils(6) = 0.0;
  permeability_stencils *= is_permeable_stencils;
  porosity_stencils *= is_permeable_stencils;
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
          SpecificHeatCapacity{capacity})};
  // logs and properties
  const Logs::Rocks::HeatLogs heat_logs{
      solid_density_stencils,
      solid_specific_heatcapacity_stencils,
      heatconductivity_stencils,
      porosity_stencils,
      water,
      grid2D->first_coord};
  const Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};
  const FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};
  // well
  const Well_KH well{
      water, core_data.is_permeable, core_data.permeability};
  // history
  const std::vector<RealType> time_steps{t1 - t0};
  const std::vector<RealType> rates{well_rate};
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

  const auto& solver{*solver_ptr};

  SolverManager solver_manager{history, solver_ptr};

  solver_manager.run(0.01);

  // assert solution
  const double tol = 1E-15;

  cout << "flow_field.axes1_as_face_normal:\n";
  cout << rates_factory.get_flow_in_axes1() << endl
       << endl;
  cout << "flow_field.axes2_as_face_normal:\n";
  cout << rates_factory.get_flow_in_axes2() << endl
       << endl;

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
    //  CHECK(v1(row, col) == v1(row + 1, col) + v2(row, col));
  }

  for (auto row{0ll}; row < grid2D->first_coord.mesh_size(); ++row)
  {
    for (auto col{0ll}; col < grid2D->second_coord.mesh_size(); ++col)
    {
      INFO("" << "col: " << col << ", row: " << row << ", bottom: " << -v1(row + 1ll, col) << ", top: " << v1(row, col) << ", right: " << v2(row, col + 1ll) << ", left: " << -v2(row, col));
      CHECK(-v1(row + 1ll, col) + v1(row, col) == v2(row, col + 1ll) - v2(row, col));
    }
  }

  const auto &[times, states] = solver.solution();
  for (size_t i{0ll}; i < times.size(); ++i)
  {
    const auto &state = states[i];
    for (auto row{0ll}; row < state.rows(); ++row)
    {
      CHECK(state(row, 0ll) <= inlet_temperature);
      for (auto col{1ll}; col < state.cols(); ++col)
      {
        CHECK(state(row, col) <= inlet_temperature);
        INFO("time: " << i << ", col: " << col << ", row: " << row);
        CHECK(state(row, col) <= state(row, col - 1ll));
      }
    }
  }

  for (size_t i{1ll}; i < times.size(); ++i)
  {
    for (auto row{0ll}; row < states[0].rows(); ++row)
    {
      for (auto col{1ll}; col < states[0].cols(); ++col)
      {
        CHECK(states[i](row, col) >= states[i - 1ull](row, col));
      }
    }
  }

  const auto precision{1e-5};
  for (auto i{times.size() - 1ll}; i < times.size(); ++i)
  {
    const auto &state = states[i];
    std::string path{std::string{"T_"} + std::to_string(0) + std::string{".txt"}};
    std::ofstream f{path};

    f << (state.cur_state / precision).round() * precision;
    f.close();
  }
}