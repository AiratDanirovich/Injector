// #include <cmath>
#include <memory>
#include <fstream>
// #include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/State2D.hpp>

#include <Injector/Grids/Factory.hpp>
#include <Injector/Solver/SolverFactory.hpp>

// #include <Injector/Properties/Logs.hpp>
// #include <Injector/Model/Phases/FluidFactory.hpp>
// #include <Injector/Model/HydrodynamicSolver.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

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

/// @brief Initial emperature is assumed to be constant
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
auto initialcondition_factory(RealType t0, const Grid_t_ptr grid, const RealType val)
{
  return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{val}, t0)};
}

struct FunctorBC : public BoundaryConditions::BCFunctorBase
{
  using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
  FunctorBC(
      RealType inlet_temp,
      const PhaseProperties &fluid,
      const Properties::ReservoirFlowField &flow_field,
      const cptr<Grid2D_t> grid_ptr)
      : inlet_temp{inlet_temp},
        fluid{fluid},
        flow_field{flow_field},
        grid_ptr{grid_ptr}
  {
  }

  RealType operator()(RealType z, RealType r, RealType t) const override
  {
    if (z == grid_ptr->first_coord.dual_front())
    {
      if (r == grid_ptr->second_coord.mesh_front())
        return fluid.volumetric_heat_capacity * flow_field.axes1_as_face_normal(0, 0) * inlet_temp;
    }

    return 0.0;
  }

  // RealType south(std::ptrdiff_t r, RealType t) const
  // {
  //   return 0.0;
  // }
  // RealType north(std::ptrdiff_t r, RealType t) const
  // {
  //   return r == 0ll ? fluid.volumetric_heat_capacity * flow_field.axes1_as_face_normal(0, 0) * inlet_temp : (RealType)0.0;
  // }

  // RealType east(std::ptrdiff_t z, RealType t) const
  // {
  //   return 0.0;
  // }
  // RealType west(std::ptrdiff_t z, RealType t) const
  // {
  //   return 0.0;
  // }

protected:
  RealType inlet_temp;
  const PhaseProperties &fluid;
  const Properties::ReservoirFlowField &flow_field;
  const cptr<Grid2D_t> grid_ptr;
};

using VR = std::vector<RealType>;

TEST_CASE("Solver", "SelfSimilarCyl")
{
  /*START*/
  // input parameters
  /*heat rate*/
  RealType q{1.0};
  /*fluid*/
  RealType viscosity{6e-4}, density{1000}, capacity{4200};
  /*collector*/
  const RealType rMin{1.0}, rMax{2.0}, zTop{0.0};
  const std::ptrdiff_t rNodes{301ull};
  const std::ptrdiff_t nLayers{11ull};
  const VR thickness(nLayers, 1); // each layer is 1m thick

  const VR conductivity(nLayers, 3.9);
  const VR is_permeable(nLayers, 1.0);
  const VR porosity(nLayers, 1e-16);
  const VR permeability(nLayers, 0.5);
  const VR solid_density(nLayers, 3.9 /*should be 2600 in SI*/);
  const VR solid_specific_heatcapacity(nLayers, 1.0 /*should be 770 in SI*/);
  /*temporal grid*/
  const std::ptrdiff_t time_steps_nmbr{501ull};
  const RealType t0{1.0}; // initial time moment
  const RealType t1{t0 + 1.0};
  const RealType time_step{(t1 - t0) / time_steps_nmbr};
  const VR time_intervals(time_steps_nmbr, time_step);
  /*temperatures*/
  const RealType well_rate{0.0};
  const RealType initial_temperature{0.0};
  const RealType inlet_temperature{1.0};
  /*END*/

  // make grid1D
  const VR z_stencils(
      Grids::Factory::generate_dual_grid_stencils_from_steps(
          zTop, thickness));
  const RealType zBottom{z_stencils.back()};
  const VR r_stencils{
      Grids::Factory::generate_dual_grid_stencils_uniform(
          Segment{rMin, rMax}, rNodes)};
  // make grid2D
  const CylinderGridFactory grid_factory{
      z_stencils, r_stencils};
  const HydrodynamicLogsFactory hydro_logs_factory{
      is_permeable, porosity, permeability,
      grid_factory.grid()->first_coord};

  // heat conductivity
  const Properties::HeatConductivity conductivity_field{
      Properties::Factory::generate_heatconductivity_Property(
          conductivity, grid_factory.grid())};
  // make fluid
  const Water water{
      FluidFactory::create_water(
          Viscosity{viscosity},
          Density{density},
          SpecificHeatCapacity{capacity})};

  SolverFactory solver_factory{
      initial_temperature,
      grid_factory};

  // volumetric heat capacity of multiphaase system
  const Properties::HeatVolumetricCapacity capacity_field{
      Properties::Factory::generate_volumetric_heatcapacity_Property(
          is_permeable, porosity,
          solid_density, solid_specific_heatcapacity,
          water, grid_factory.grid())};

  const Well_KH well{
      Phases::FluidFactory::create_water(0.0, 0.0),
      hydro_logs_factory.is_permeable,
      hydro_logs_factory.permeability};
  Properties::ReservoirFlowField flow_field{
      well_rate, well, *grid_factory.grid()};

  // initial conditions
  const auto initial_state{initialcondition_factory(t0, grid_factory.grid(), initial_temperature)};
  // boundary conditions
  const GPN::BoundaryConditions::BoundaryConditions bc{
      *grid_factory.grid(),
      std::make_shared<FunctorBC>(
          inlet_temperature, water, flow_field, grid_factory.grid()),
      BoundaryConditions::BoundaryCondition::second};
  // time moments
  const VR t_stencils(
      Grids::Factory::generate_dual_grid_stencils_from_steps(
          t0, time_intervals));
  // solver

  Solver solver{
      conductivity_field,
      flow_field, grid_factory.grid(),
      solver_factory.capacity_field(),
      initial_state,
      bc, t0};

  const double tol = 1E-3;

  // assert initial condition
  // for (std::ptrdiff_t col{0ll}; col < initial_state.cols(); ++col)
  // {
  //   for (std::ptrdiff_t row{0ll}; row < initial_state.rows(); ++row)
  //   {
  //     const auto [z, r] = grid2D->coordinates(row, col);
  //     const auto val{es(z, r, t0)};
  //     INFO("" << "col: " << col << ", row: " << row << ", calc: " << initial_state(row, col) << ", ref: " << val);
  //     CHECK_THAT(initial_state(row, col), WithinRel(val, tol));
  //   }
  // }

  std::string pathr{"data_r.csv"};
  std::string pathz{"data_z.csv"};

  // assert solution
  for (size_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
  {
    solver.advance(time_intervals[t_step]);
  }

  std::ofstream fr{pathr};
  assert(fr.is_open());
  fr << "r;Tcalc\n";
  const auto &[time, state] = solver.solution().back();
  for (std::ptrdiff_t col{0ll}; col < state.cols(); ++col)
  {
    for (std::ptrdiff_t row{nLayers / 2}; row < nLayers / 2 + 1ll /*state.rows()*/; ++row)
    {
      const auto [z, r] = grid_factory.grid()->coordinates(row, col);

      fr << r << ';' << state(row, col) << '\n';
    }
  }
  fr.close();

  // std::ofstream fz{pathz};
  // assert(fz.is_open());
  // fz << "z;TcalcL;TcalcM;TcalcR;Tref\n";
  // for (std::ptrdiff_t row{0}; row < state.rows(); ++row)
  // {
  //   for (std::ptrdiff_t col{rNodes/2}; col < rNodes/2+1ll /*state.cols()*/; ++col)
  //   {
  //     const auto [z, r] = grid2D->coordinates(row, col);
  //     const auto val{es(z, r, time)};
  //     const auto val2{es(row, col, time, *grid2D)};

  //     assert(val == val2);

  //     fz
  //     << z << ';'
  //     << val << ';'
  //     << state(row, 0)<< ';'
  //     << state(row, col)<< ';'
  //     << state(row, state.cols()-1ll)
  //     << '\n';

  //     INFO("" << "col: " << col << ", row: " << row << ", calc: " << state(row, col) << ", ref: " << val);
  //     CHECK_THAT(state(row, col), WithinRel(val, tol));
  //   }
  // }
  // fz.close();
}