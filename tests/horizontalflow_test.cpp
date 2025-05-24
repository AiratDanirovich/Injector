// #include <cmath>
#include <memory>
#include <fstream>
#include <string>

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/State2D.hpp>

#include <Injector/Grids/Factory.hpp>
#include <Injector/Solver/SolverFactory.hpp>

// #include <Injector/Properties/Logs.hpp>
// #include <Injector/Model/Phases/FluidFactory.hpp>
// #include <Injector/Model/HydrodynamicSolver.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

#include <Injector/Properties/FieldFactory.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::Logs;
using namespace GPN::Grids;
using namespace GPN::Phases;
// using namespace GPN::Model::Injector;
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
      const PhaseProperties &fluid, // inlet fluid
      const Properties::ReservoirFlowField &flow_field, // volumetric flow rate
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
  const Properties::ReservoirFlowField &flow_field;
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
const VR is_permeable(nLayers, 1.0);
const VR porosity(nLayers, 1.0);
const VR permeability(nLayers, 0.5);
const VR ext_pressure(nLayers, 13E6);
const VR skin(nLayers, 0.0);

// heat logs
const VR conductivity(nLayers, 3.9);
const VR solid_density(nLayers, 3.9 /*should be 2600 in SI*/);
const VR solid_specific_heatcapacity(nLayers, 1.0 /*should be 770 in SI*/);
/*temporal grid*/
const std::ptrdiff_t time_steps_nmbr{5ull};
const RealType t0{1.0}; // initial time moment
const RealType t1{t0 + 0.01};
const RealType time_step{(t1 - t0) / time_steps_nmbr};
const VR time_intervals(time_steps_nmbr, time_step);
/*temperatures*/
const RealType well_rate{0.01};
const RealType initial_temperature{0.0};
const RealType inlet_temperature{1.0};
/*END*/

TEST_CASE("Solver", "SelfSimilarCyl")
{
  // make grid1D
  const VR z_stencils(
      Grids::Factory::generate_dual_grid_stencils_from_steps(
          zTop, thickness));
  const RealType zBottom{z_stencils.back()};
  const VR r_stencils{
      Grids::Factory::generate_dual_grid_stencils_uniform(
          Segment{rMin, rMax}, rNodes)};
  // make grid2D
  const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
  const HydrodynamicLogsFactory hydro_logs_factory{
      is_permeable, porosity, permeability, ext_pressure, skin,
      grid2D->first_coord};

  // heat conductivity
  const Properties::HeatConductivity conductivity_field{
      Properties::Factory::generate_heatconductivity_Property(
          conductivity, grid2D)};

  // make fluid
  const Water water{
      FluidFactory::create_water(
          Viscosity{viscosity},
          Density{density},
          SpecificHeatCapacity{capacity})};

  SolverFactory solver_factory{
      initial_temperature,
      grid2D};

  // volumetric heat capacity of multiphase system
  const auto capacity_field{
      Properties::Factory::generate_volumetric_heatcapacity_Property(
          is_permeable, porosity,
          solid_density, solid_specific_heatcapacity,
          water, grid2D)};

  // // const Well_KH well{
  // //     Phases::FluidFactory::create_water(0.0, 0.0),
  // //     hydro_logs_factory.is_permeable,
  // //     hydro_logs_factory.permeability};
  // Properties::ReservoirFlowField flow_field{
  //   Properties::FlowFactory::horizontal_flow(well_rate*capacity_field.value(0,0), hydro_logs_factory.is_permeable, *grid2D)
  //   //  well_rate, well, *grid2D
  //   };

  // // initial condition
  // const auto initial_state{ICFactory(t0, grid2D, initial_temperature)};
  // // boundary conditions
  // const GPN::BoundaryConditions::BoundaryConditions bc{
  //     *grid2D,
  //     std::make_shared<FunctorBC>(
  //         inlet_temperature, water, flow_field, grid2D),
  //     BoundaryConditions::BoundaryCondition::second};
  // // time moments
  // const VR t_stencils(
  //     Grids::Factory::generate_dual_grid_stencils_from_steps(
  //         t0, time_intervals));
  // // solver

  // Solver solver{
  //     conductivity_field,
  //     flow_field, grid2D,
  //     solver_factory.capacity_field(),
  //     initial_state,
  //     bc, t0};

  // std::string pathr{"data_r.csv"};
  // std::string pathz{"data_z.csv"};

  // // assert solution
  // const double tol = 1E-15;
  // for (size_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
  // {
  //   // for (auto i{0ull}; i < matricies.mx.size(); ++i)
  //   // {
  //   //   const auto &m = matricies.mx[i];
  //   //   REQUIRE(m.rows() == m.cols());
  //   //   {
  //   //     auto col{0ll};
  //   //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << 0.0 << ", r: " << m.coeff(col, col + 1ll));
  //   //     CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col + 1ll), tol));
  //   //   }
  //   //   for (auto col{1ll}; col < m.cols() - 1; ++col)
  //   //   {
  //   //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << m.coeff(col, col - 1ll) << ", r: " << m.coeff(col, col + 1ll));
  //   //     CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll) + m.coeff(col, col + 1ll), tol));
  //   //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col));
  //   //     CHECK(m.coeff(col, col) > tol);
  //   //   }
  //   //   {
  //   //     auto col{m.cols() - 1ll};
  //   //     INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << m.coeff(col, col - 1ll) << ", r: " << 0.0);
  //   //     CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll), tol));
  //   //   }
  //   // }
  // }

  // const auto &[times, states] = solver.solution();

  // for (auto i{0ull}; i < times.size(); ++i)
  // {
  //   std::string path{std::string{"T_"} + std::to_string(i) + std::string{".txt"}};
  //   std::ofstream f{path};

  //   f << states[i].cur_state;
  //   f.close();
  // }
}