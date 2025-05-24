#include <iostream>

#include <Injector/Grids/Defines.h>
#include <Injector/Solver/SolverFactory.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

TEST_CASE("Solver")
{
  const double tol = 4E-11;

  RealType val{1.0};
  RealType well_rate{0.0};
  SolverFactory solver_factory{val, Box{Segment{0, 1}, Segment{0, 1}}, 31, 51};

  const auto grid2D{solver_factory.grid()};
  const auto conductivity_field{solver_factory.conductivity_field()};
  const auto capacity_field{solver_factory.capacity_field()};
  const auto initial_state{solver_factory.initial_state()};
  const auto bc{solver_factory.boundary_conditions()};

  const auto flow_field{solver_factory.flow_field(well_rate)};

  Solver solver{
      conductivity_field,
      flow_field, grid2D,
      capacity_field,
      initial_state,
      bc, 0.0};

  RealType step{1.0};
  ptrdiff_t nT{3};
  for (ptrdiff_t t{0ll}; t < nT; ++t)
  {
    solver.advance(step);
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
