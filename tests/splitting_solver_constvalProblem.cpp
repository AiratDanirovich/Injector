#include <iostream>
#include <memory>

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/History/ZeroRatesFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace GPN;

struct ABCFunctor : public GPN::BoundaryConditions::BCFunctorBase
{
  ABCFunctor(RealType val) : val{val} {}

  RealType operator()(const ptrdiff_t z, const RealType r, const RealType t) const override
  {

    return val;
  }

  RealType operator()(const RealType z, const ptrdiff_t r, const RealType t) const override
  {
    return val;
  }

protected:
  RealType val;
};

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

using VR = std::vector<RealType>;

const auto box{Box{Segment{0, 1}, Segment{0, 1}}};
const auto nLayers{31ll}, nR{51ll};

const VR is_permeable_stencils(nLayers - 1ll, 1.0);
const LogValuesContainer porosity_stencils{LogValuesContainer::Constant(nLayers - 1ll, 1e-16)};
const LogValuesContainer ext_pressure_stencils{LogValuesContainer::Constant(nLayers-1ll, 260 * 1e5)};

const LogValuesContainer solid_heatconductivity_stencils(LogValuesContainer::Constant(nLayers - 1ll, 3.9));
const LogValuesContainer solid_density_stencils{LogValuesContainer::Constant(nLayers - 1ll, 3.9 /*should be 2600 in SI*/)};
const LogValuesContainer solid_specific_heatcapacity_stencils{LogValuesContainer::Constant(nLayers - 1ll, 1.0 /*should be 770 in SI*/)};

TEST_CASE("Solver")
{
  const double tol = 4E-11;

  RealType val{1.0};
  RealType well_rate{0.0};

  const auto grid2D{Grids::CylinderGridFactory::create(box, nLayers, nR)};
  const auto &grid_z{grid2D->first_coord()};

  const State::State2D initial_state{State::State2D::FillWithConst(*grid2D, val)};

  const BoundaryConditions::BoundaryConditions bc{
      *grid2D,
      std::make_shared<ABCFunctor>(ABCFunctor{val})};

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
      Phases::FluidFactory::create_water(1.0, 1.0),
      grid_z};

  const Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  const FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};

  // external pressure log
  const auto external_pressure{
      Logs::ExtPressureFactory::create(
          ext_pressure_stencils,
          is_permeable_stencils,
          grid_z)};
  // rates field factory
  FaceProperties::ZeroRatesFactory rates_factory{
      grid2D, core_data.is_permeable, external_pressure};

  Solver solver{
      heat_face_props.medium_heat_conductivity,
      grid2D,
      heat_props.medium_vol_heatcapacity,
      rates_factory, initial_state,
      bc, 0.0};

  const RealType step{1.0};
  const ptrdiff_t nT{3};
  for (ptrdiff_t t{0ll}; t < nT; ++t)
  {
    solver.advance(step);
    solver.save_state();
    std::cout << "time: " << t * step << std::endl;
  }

  auto [time, solution] = solver.solution().back();

  CHECK(time == nT * step);
  for (std::ptrdiff_t col{0ll}; col < solution.cols(); ++col)
    for (std::ptrdiff_t row{0ll}; row < solution.rows(); ++row)
    {
      INFO("" << "row: " << row << ", col: " << col << ", solution: " << solution(row, col));
      CHECK_THAT(solution(row, col), WithinRel(val, tol));
    }
}
