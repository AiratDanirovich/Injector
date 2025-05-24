#include <cmath>
#include <memory>
#include <fstream>
#include <cassert>

#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/HydrodynamicSolver.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>
#include <Injector/Solver/SplittingMethod/SolverFactory.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::Logs;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::Model::Injector;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

struct ExactSolution
{
  ExactSolution(
      const Properties::HeatVolumetricCapacity &volumetric_capacity,
      const Properties::HeatConductivity &heat_conductivity,
      RealType q)
      : kappa{volumetric_capacity, heat_conductivity},
        heat_conductivity{heat_conductivity},
        grid{*heat_conductivity.grid},
        q{q}
  {
  }

  RealType operator()(RealType z, RealType r, RealType t) const
  {
    return -q * heat_conductivity.value(0ull, 0ull) *
           std::expint(-r * r / (4.0 * kappa.value(0ull, 0ull) * t));
  }

  template <typename Grid_t>
  RealType operator()(ptrdiff_t z, ptrdiff_t r, RealType t, const Grid_t &grid) const
  {
    const auto [zv, rv]{grid.coordinates(z, r)};
    return (*this)(zv, rv, t);
  }

protected:
  const Properties::ThermalDiffusivity kappa;
  const Properties::HeatConductivity &heat_conductivity;
  const RealType q;
  const Properties::HeatConductivity::Grid_type &grid;
};

struct FunctorIC
{
  FunctorIC(const ExactSolution &es)
      : es{es}
  {
  }
  RealType operator()(RealType z, RealType r, RealType t0) const
  {
    return es(z, r, t0);
  }

protected:
  const ExactSolution &es;
};

template <typename Grid_t_ptr>
auto initialcondition_factory(RealType t0, const Grid_t_ptr grid, const ExactSolution &es)
{
  return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{es}, t0)};
}

struct FunctorBC : public GPN::BoundaryConditions::BCFunctorBase
{
  FunctorBC(const ExactSolution &es)
      : es{es}
  {
  }
  RealType operator()(RealType z, RealType r, RealType t) const override
  {
    return es(z, r, t);
  }

protected:
  const ExactSolution &es;
};

using VR = std::vector<RealType>;

TEST_CASE("SolverManager", "SelfSimilarCyl")
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
  const VR thickness(nLayers, 0.01); // each layer is 1m thick

  const VR conductivity(nLayers, 3.9);
  const VR porosity(nLayers, 1e-16);
  const VR is_permeable(nLayers, 1.0);
  const VR solid_density(nLayers, 3.9 /*should be 2600 in SI*/);
  const VR solid_specific_heatcapacity(nLayers, 1.0 /*should be 770 in SI*/);
  /*temporal grid*/
  const std::ptrdiff_t time_steps_nmbr{501ull};
  const RealType t0{1.0}; // initial time moment
  const RealType t1{t0 + 1.0};
  const auto time_moments{
      Grids::Factory::generate_dual_grid_stencils_uniform(
          t0, t1, time_steps_nmbr + 1)};
  const RealType time_step{(t1 - t0) / time_steps_nmbr};
  // const VR time_intervals(time_steps_nmbr, time_step);
  const VR time_intervals{
      Grids::Factory::generate_dual_grid_steps(
          time_moments)};
  /*rates*/
  const VR rates{Logs::Factory::generate_rates_StepProperty(
      time_moments)};
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
  const auto grid2D{
      Grids::Factory::create_cylinder_grid_2D_ptr(
          z_stencils, r_stencils)};
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
  // volumetric heat capacity of multiphaase system
  const Properties::HeatVolumetricCapacity capacity_field{
      Properties::Factory::generate_volumetric_heatcapacity_Property(
          is_permeable, porosity,
          solid_density, solid_specific_heatcapacity,
          water, grid2D)};
  // exact solution
  ExactSolution es{capacity_field,
                   conductivity_field, q};
  // initial conditions
  const auto initial_state{initialcondition_factory(t0, grid2D, es)};
  // boundary conditions
  const GPN::BoundaryConditions::BoundaryConditions bc{
      *grid2D, std::make_shared<FunctorBC>(es)};
  // history
  const History history{
      InjectorRate{
          StepPropertyGrid{
              rates,
              TemporalGridDual{GridDualStencils{DualStepsContainer{time_intervals}}}}}};
  // const VR t_stencils(
  //     Grids::Factory::generate_dual_grid_stencils_from_steps(
  //         t0, time_intervals));

  // solver
  using Solver_t = Solver<
      Grids::StructuredCylinderGrid2DAxisymmetric,
      Properties::HeatVolumetricCapacity,
      FlowField>;

  FlowField flow_field{};

  const cptr<Solver_t> solver{std::make_shared<Solver_t>(
      conductivity_field,
      flow_field, grid2D,
      capacity_field,
      initial_state,
      bc, t0)};

  // const double tol = 1E-3;

  SolverManager manager{
      history,
      solver};

  manager.run(time_intervals[0] / 2.5);
}