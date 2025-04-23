#include <cmath>
#include <memory>
#include <fstream>
#include <cassert>

#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>
#include <Injector/Solver/SplittingMethod/SolverFactory.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::Phases;
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
  
  template<typename Grid_t>
  RealType operator()(ptrdiff_t z, ptrdiff_t r, RealType t, const Grid_t& grid) const
  {
    const auto[zv, rv]{grid.coordinates(z,r)};
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

TEST_CASE("Solver", "SelfSimilarCyl")
{
  /*START*/
  // input parameters
  /*heat rate*/
  RealType q{1.0};
  /*fluid*/
  RealType viscosity{6e-4}, density{1000}, capacity{4200};
  /*collector*/
  const RealType rMin{2.0}, rMax{3.0}, zTop{0.0};
  std::size_t rNodes{101ull};
  const std::size_t nLayers{21ull};
  const VR thickness(nLayers, 0.01); // each layer is 1m thick

  const VR conductivity(nLayers, 3.9);
  const VR porosity(nLayers, 1e-16);
  const VR is_permeable(nLayers, 1.0);
  const VR solid_density(nLayers, 3.9 /*should be 2600 in SI*/);
  const VR solid_specific_heatcapacity(nLayers, 1.0 /*should be 770 in SI*/);
  /*temporal grid*/
  const std::size_t time_steps_nmbr{11ull};
  const RealType t0{1.0}; // initial time moment
  const RealType t1{t0 + 1.0};
  const RealType time_step{(t1 - t0) / time_steps_nmbr};
  const VR time_intervals(time_steps_nmbr, time_step);
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
  // time moments
  const VR t_stencils(
      Grids::Factory::generate_dual_grid_stencils_from_steps(
          t0, time_intervals));
  // solver
  Solver solver{
      conductivity_field, grid2D,
      capacity_field,
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

  std::string path{"data.csv"};

  std::ofstream f{path};
  assert(f.is_open());
  f << "r;Tref;Tcalc\n";

  // assert solution
  for (size_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
  {
    solver.advance(time_intervals[t_step]);
  }

  const auto &[time, state] = solver.solution().back();
  for (std::ptrdiff_t col{0ll}; col < state.cols(); ++col)
  {
    for (std::ptrdiff_t row{nLayers/2}; row < nLayers/2+1ll /*state.rows()*/; ++row)
    {
      const auto [z, r] = grid2D->coordinates(row, col);
      const auto val{es(z, r, time)};
      const auto val2{es(row, col, time, *grid2D)};

      assert(val == val2);
      
      f << r << ';' << val << ';' << state(row, col) << '\n';

      INFO("" << "col: " << col << ", row: " << row << ", calc: " << state(row, col) << ", ref: " << val);
      CHECK_THAT(state(row, col), WithinRel(val, tol));
    }
  }
}