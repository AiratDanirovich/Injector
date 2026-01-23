#include <cmath>
#include <memory>
#include <fstream>
#include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/Factory.hpp>
#include <Injector/History/History.hpp>
#include <Injector/History/ZeroRatesFactory.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/generate_stencils_and_steps.hpp"

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

struct ExactSolution
{
  using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;

  ExactSolution(
      const Properties::MediumHeatVolumetricCapacity<Grid2D_t> &volumetric_capacity,
      const Properties::HeatConductivity<Grid2D_t> &heat_conductivity,
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

  const Grid2D_t &grid;

protected:
  const Properties::ThermalDiffusivity<Grid2D_t> kappa;
  const Properties::HeatConductivity<Grid2D_t> &heat_conductivity;
  const RealType q;
};

struct FunctorIC : public InitialConditions::ICFunctorBase
{
  FunctorIC(const ExactSolution &es)
      : es{es}
  {
  }
  RealType operator()(const ptrdiff_t z_id, const ptrdiff_t r_id, const RealType t0) const override
  {
    const RealType z = es.grid.first_coord().mesh_nodes(z_id);
    const RealType r = es.grid.second_coord().mesh_nodes(r_id);
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

struct FunctorBC : public GPN::BoundaryConditions::GeneralBC::BCFunctorBase
{
  using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;

  FunctorBC(
      const ExactSolution &es,
      const cptr<Grid2D_t> grid2D)
      : es{es},
        grid2D{grid2D}
  {
  }
  BC_descriptor operator()(const ptrdiff_t z_id, RealType r, const RealType t,
                            const BCType) const override
  {
    RealType z{grid2D->first_coord().mesh_nodes(z_id)};
    if (r == grid2D->second_coord().dual_front())
      r = grid2D->second_coord().mesh_front();
    else if (r == grid2D->second_coord().dual_back())
      r = grid2D->second_coord().mesh_back();
    else
      assert(false);
    return BC_descriptor::BC_I(es(z, r, t));
  //  return es(z, r, t);
  }

  BC_descriptor operator()(RealType z, const ptrdiff_t r_id, const RealType t,
                            const BCType) const override
  {
    RealType r{grid2D->second_coord().mesh_nodes(r_id)};
    if (z == grid2D->first_coord().dual_front())
      z = grid2D->first_coord().mesh_front();
    else if (z == grid2D->first_coord().dual_back())
      z = grid2D->first_coord().mesh_back();
    else
      assert(false);
    return BC_descriptor::BC_I(es(z, r, t));
  //  return es(z, r, t);
  }

protected:
  const ExactSolution &es;
  const cptr<Grid2D_t> grid2D;
};

using VR = std::vector<RealType>;

TEST_CASE("SolverManager", "SelfSimilarCyl")
{
  /*START*/
  // input parameters
  /*heat rate*/
  RealType q{1.0};
  /*fluid*/
  RealType viscosity{6e-4}, density{1000}, capacity{4200}, heat_conductivity{0.6};
  /*collector*/
  const RealType rMin{1.0}, rMax{2.0}, zTop{0.0};
  const std::ptrdiff_t rNodes{101ull};
  const std::ptrdiff_t nLayers{3ull};
  const VR thickness(nLayers, 0.01); // each layer is 1m thick

  const VR porosity(nLayers, 1e-16);
  const VR is_permeable_stencils(nLayers, 1.0);
  const LogValuesContainer porosity_stencils{LogValuesContainer::Constant(nLayers, 0.0)};
  const LogValuesContainer ext_pressure_stencils{LogValuesContainer::Constant(nLayers, 260 * 1e5)};
  const LogValuesContainer solid_heatconductivity_stencils{LogValuesContainer::Constant(nLayers, 3.9)};
  const LogValuesContainer solid_density_stencils{LogValuesContainer::Constant(nLayers, 3.9 /*should be 2600 in SI*/)};
  const LogValuesContainer solid_specific_heatcapacity_stencils{LogValuesContainer::Constant(nLayers, 1.0 /*should be 770 in SI*/)};
  /*temporal grid*/
  const std::ptrdiff_t time_steps_nmbr{15ull};
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
  // const VR rates{Logs::RawDataFactory::generate_rates(
  //     time_moments)};
  /*END*/

  // make grid1D
  const auto z_stencils(
      Grids::Factory::generate_dual_grid_stencils_from_steps(
          zTop, thickness));
  const auto r_stencils{
      Grids::Factory::generate_dual_grid_stencils_uniform(
          Segment{rMin, rMax}, rNodes)};
  const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
  const auto &grid_z{grid2D->first_coord()};
  // make fluid
  const PhaseProperties water{
      FluidFactory::create_water(
          Viscosity{viscosity},
          Density{density},
          SpecificHeatCapacity{capacity},
          GPN::HeatConductivity{heat_conductivity})};
  // solver
  const Logs::Rocks::IsPermeableLog core_data{
      is_permeable_stencils,
      grid_z};

  const auto flow_field{
      FaceProperties::FlowFactory::zero_flow(
          core_data.is_permeable, *grid2D)};

  const Logs::Rocks::HeatLogs heat_logs{
      solid_density_stencils,
      solid_specific_heatcapacity_stencils,
      solid_heatconductivity_stencils,
      porosity_stencils,
      water,
      grid_z};

  const Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  const FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};
  // history
  const std::vector<RealType> rates(time_intervals.size(), 1.0);
  const std::vector<RealType> inlet_temperature_set(time_intervals.size(), 0.0);
  const ptr<History> history{make_shared<History>(
      HistoryFactory::createFixedRate(time_intervals, rates, inlet_temperature_set))};

  // exact solution
  ExactSolution es{heat_props.medium_vol_heatcapacity,
                   heat_props.medium_heat_conductivity_axes1, q};
  // initial conditions
  const auto initial_state{initialcondition_factory(t0, grid2D, es)};
  // boundary conditions
  const GPN::BoundaryConditions::GeneralBC bc{
      grid2D, std::make_shared<FunctorBC>(es, grid2D),
        BoundaryConditions::GeneralBC::BoundaryCondition::first};
  // external pressure log
  const auto external_pressure{
      Logs::ExtPressureFactory::create(
          ext_pressure_stencils,
          is_permeable_stencils,
          grid_z)};
  // rates field factory
  FaceProperties::ZeroRatesFactory rates_factory{
      grid2D, core_data.is_permeable, external_pressure};

  // solver
  using Solver_t = decltype(Solver{
      heat_face_props.medium_heat_conductivity,
      grid2D,
      heat_props.medium_vol_heatcapacity,
      rates_factory, initial_state,
      bc, t0});

  const ptr<Solver_t> solver{std::make_shared<Solver_t>(
      heat_face_props.medium_heat_conductivity,
      grid2D,
      heat_props.medium_vol_heatcapacity,
      rates_factory, initial_state,
      bc, t0)};

  SolverManager manager{
      history,
      solver};

  manager.run(time_intervals[0] / 2.5);
}