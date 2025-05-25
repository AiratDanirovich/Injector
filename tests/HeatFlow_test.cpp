#include <memory>
#include <fstream>
#include <string>
#include <numbers>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

#include <Injector/Properties/FieldsFactory.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

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
  FunctorBC(
      RealType inlet_temp,
      const Logs::IsPermeable &is_permeable,
      const FaceProperties::ReservoirFlowField &flow_field, // volumetric flow rate
      const cptr<const Grid2D_t> grid_ptr)
      : inlet_temp{inlet_temp},
        flow_field{flow_field},
        is_permeable{is_permeable},
        grid_ptr{grid_ptr}
  {
  }

  RealType operator()(ptrdiff_t z_id, RealType r, RealType t) const override
  {
    return 0.0;
  }

  RealType operator()(RealType z, ptrdiff_t r_id, RealType t) const override
  {
    if (z == grid_ptr->first_coord.dual_front())
      return flow_field.axes1_as_face_normal.topRows(1ll)(1ll, r_id) * inlet_temp;

    if (z == grid_ptr->first_coord.dual_back())
      return flow_field.axes1_as_face_normal.bottomRows(1ll)(1ll, r_id) * inlet_temp;

    return 0.0;
  }

protected:
  RealType inlet_temp;
  const FaceProperties::ReservoirFlowField &flow_field;
  const Logs::IsPermeable &is_permeable;
  const cptr<const Grid2D_t> grid_ptr;
};

using VR = std::vector<RealType>;

/*START*/
// input parameters
/*fluid*/
RealType viscosity{6e-4}, density{1}, capacity{1};
/*collector*/
const RealType rMin{1 / (2 * numbers::pi)}, rMax{1.0}, zTop{0.0}; // m
const std::ptrdiff_t rNodes{201ull};
const std::ptrdiff_t nLayers{5ull};
const VR thickness(nLayers, 1.0); // each layer is 1m thick

// hydrodynamic logs
VR is_permeable_stencils(nLayers, 1.0);
LogValuesContainer porosity_stencils{LogValuesContainer::Constant(nLayers, 1.0)};
const VR permeability_stencils(nLayers, 0.5);

// heat logs
const VR heatconductivity_stencils(nLayers, 0.0);
const LogValuesContainer solid_density_stencils{LogValuesContainer::Constant(nLayers, 3.9 /*should be 2600 in SI*/)};
const LogValuesContainer solid_specific_heatcapacity_stencils{LogValuesContainer::Constant(nLayers, 1.0 /*should be 770 in SI*/)};
/*temporal grid*/
const std::ptrdiff_t time_steps_nmbr{51ull};
const RealType t0{0.0};      // initial time moment
const RealType t1{t0 + 0.1}; // s
const RealType time_step{(t1 - t0) / time_steps_nmbr};
const VR time_intervals(time_steps_nmbr, time_step);
const VR t_stencils(
    Grids::Factory::generate_dual_grid_stencils_from_steps(
        t0, time_intervals));
/*temperatures*/
const RealType well_rate{1}; // m^3/s
const RealType initial_temperature{0.0};
const RealType inlet_temperature{1.0};
/*END*/

TEST_CASE("Solver", "SelfSimilarCyl")
{
  // make grid2D
  const auto grid2D{
      Grids::CylinderGridFactory::create(
          Grids::Factory::generate_dual_grid_stencils_from_steps(
              zTop, thickness),
          Grids::Factory::generate_dual_grid_stencils_uniform(
              Segment{rMin, rMax}, rNodes))};
  const auto &grid{grid2D->first_coord};

  const Logs::IsPermeable is_permeable{
      IsPermeableFactory::create(is_permeable_stencils, grid)};
  const Logs::Porosity porosity{
      PorosityFactory::create(porosity_stencils, is_permeable_stencils, grid)};

  // make fluid
  const Water water{
      FluidFactory::create_water(
          Viscosity{viscosity},
          Density{density},
          SpecificHeatCapacity{capacity})};

  const Logs::Rocks::IsPermeableLog hydrodynamics_logs{
      is_permeable_stencils,
      grid2D->first_coord};

  const Logs::Rocks::HeatLogs heat_logs{
      solid_density_stencils,
      solid_specific_heatcapacity_stencils,
      heatconductivity_stencils,
      porosity.log_vals,
      water,
      grid2D->first_coord};

  const Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  const FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};

  FaceProperties::ReservoirFlowField flow_field{
      FaceProperties::FlowFactory::horizontal_flow(
          well_rate, hydrodynamics_logs.is_permeable, *grid2D)};

  FaceProperties::multiply(flow_field, water.volumetric_heat_capacity);

  // initial condition
  const auto initial_state{ICFactory(t0, grid2D, initial_temperature)};
  // boundary conditions
  const GPN::BoundaryConditions::BoundaryConditions bc{
      *grid2D,
      std::make_shared<FunctorBC>(
          inlet_temperature, is_permeable, flow_field, grid2D),
      BoundaryConditions::BoundaryCondition::second};
  // solver

  Solver solver{
      heat_face_props.heat_conductivity,
      flow_field, grid2D,
      heat_props.medium_vol_heatcapacity,
      initial_state,
      bc, t0};

  // assert solution
  const double tol = 1E-15;
  for (size_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
  {
    solver.advance(time_intervals[t_step]);
  }

  const auto &v1 = flow_field.axes1_as_face_normal;
  for (auto col{0ll}; col < v1.cols(); ++col)
    for (auto row{0ll}; row < v1.rows(); ++row)
      CHECK(v1(row, col) == 0.0);

  const auto &v2 = flow_field.axes2_as_face_normal;
  for (auto col{0ll}; col < v2.cols(); ++col)
    for (auto row{0ll}; row < v2.rows(); ++row)
      CHECK(((v2(row, col) > 0.0) || ((v2(row, col) == 0.0) && (is_permeable.log_vals(row) == 0.0))));

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
        CHECK(state(row, col) <= state(row, col - 1ll));
      }
    }
  }

  for (auto i{0ull}; i < times.size(); ++i)
  {
    std::string path{std::string{"T_"} + std::to_string(i) + std::string{".txt"}};
    std::ofstream f{path};

    f << states[i].cur_state;
    f.close();
  }
}