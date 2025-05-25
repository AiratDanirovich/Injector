#include <memory>
#include <fstream>
#include <string>

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
      const PhaseProperties &fluid,                         // inlet fluid
      const FaceProperties::ReservoirFlowField &flow_field, // volumetric flow rate
      const cptr<const Grid2D_t> grid_ptr)
      : inlet_temp{inlet_temp},
        fluid{fluid},
        flow_field{flow_field},
        grid_ptr{grid_ptr}
  {
  }

  RealType operator()(RealType z, RealType r, RealType t) const override
  {
    if (r == grid_ptr->second_coord.dual_front())
      return fluid.volumetric_heat_capacity * flow_field.axes2_as_face_normal(0, 0) * inlet_temp;

    return 0.0;
  }

protected:
  RealType inlet_temp;
  const PhaseProperties &fluid;
  const FaceProperties::ReservoirFlowField &flow_field;
  const cptr<const Grid2D_t> grid_ptr;
};

using VR = std::vector<RealType>;

/*START*/
// input parameters
/*fluid*/
RealType viscosity{6e-4}, density{1000}, capacity{4200};
/*collector*/
const RealType rMin{100.0}, rMax{101.0}, zTop{0.0};
const std::ptrdiff_t rNodes{31ull};
const std::ptrdiff_t nLayers{11ull};
const VR thickness(nLayers, 1); // each layer is 1m thick

// hydrodynamic logs
const VR is_permeable_stencils(nLayers, 1.0);
const LogValuesContainer porosity_stencils{LogValuesContainer::Constant(nLayers, 1.0)};
const VR permeability_stencils(nLayers, 0.5);
const VR ext_pressure_stencils(nLayers, 13E6);
const VR skin_stencils(nLayers, 0.0);

// heat logs
const VR heatconductivity_stencils(nLayers, 3.9);
const LogValuesContainer solid_density_stencils{LogValuesContainer::Constant(nLayers,3.9 /*should be 2600 in SI*/)};
const LogValuesContainer solid_specific_heatcapacity_stencils{LogValuesContainer::Constant(nLayers, 1.0 /*should be 770 in SI*/)};
/*temporal grid*/
const std::ptrdiff_t time_steps_nmbr{5ull};
const RealType t0{1.0}; // initial time moment
const RealType t1{t0 + 0.01};
const RealType time_step{(t1 - t0) / time_steps_nmbr};
const VR time_intervals(time_steps_nmbr, time_step);
const VR t_stencils(
    Grids::Factory::generate_dual_grid_stencils_from_steps(
        t0, time_intervals));
/*temperatures*/
const RealType well_rate{0.01};
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

  const Logs::Rocks::CoreSampleLogs core_data{
      is_permeable_stencils,
      porosity_stencils,
      permeability_stencils,
      grid2D->first_coord};

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
      porosity_stencils,
      water,
      grid2D->first_coord};

  const Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  const FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};

  FaceProperties::ReservoirFlowField flow_field{
      FaceProperties::FlowFactory::horizontal_flow(
          well_rate, hydrodynamics_logs.is_permeable, *grid2D)};

  // initial condition
  const auto initial_state{ICFactory(t0, grid2D, initial_temperature)};
  // boundary conditions
  const GPN::BoundaryConditions::BoundaryConditions bc{
      *grid2D,
      std::make_shared<FunctorBC>(
          inlet_temperature, water, flow_field, grid2D),
      BoundaryConditions::BoundaryCondition::second};
  // solver

  Solver solver{
      heat_face_props.heat_conductivity,
      flow_field, grid2D,
      heat_props.medium_vol_heatcapacity,
      initial_state,
      bc, t0};

  std::string pathr{"data_r.csv"};
  std::string pathz{"data_z.csv"};

  // assert solution
  const double tol = 1E-15;
  for (size_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
  {
    // for (auto i{0ull}; i < matricies.mx.size(); ++i)
    // {
    //   const auto &m = matricies.mx[i];
    //   REQUIRE(m.rows() == m.cols());
    //   {
    //     auto col{0ll};
    //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << 0.0 << ", r: " << m.coeff(col, col + 1ll));
    //     CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col + 1ll), tol));
    //   }
    //   for (auto col{1ll}; col < m.cols() - 1; ++col)
    //   {
    //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << m.coeff(col, col - 1ll) << ", r: " << m.coeff(col, col + 1ll));
    //     CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll) + m.coeff(col, col + 1ll), tol));
    //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col));
    //     CHECK(m.coeff(col, col) > tol);
    //   }
    //   {
    //     auto col{m.cols() - 1ll};
    //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << m.coeff(col, col - 1ll) << ", r: " << 0.0);
    //     CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll), tol));
    //   }
    // }
  }

  const auto &[times, states] = solver.solution();

  for (auto i{0ull}; i < times.size(); ++i)
  {
    std::string path{std::string{"T_"} + std::to_string(i) + std::string{".txt"}};
    std::ofstream f{path};

    f << states[i].cur_state;
    f.close();
  }
}