#include <iostream>

#include <Injector/Solver/SplittingMethod/SolverFactory.hpp>

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

  Solver solver{Factory::make_solver(val)};

  RealType step{1.0};
  ptrdiff_t nT{3};
  for (ptrdiff_t t{0ll}; t < nT; ++t)
  {
    solver.advance(step);
    std::cout << "time: " << t*step << std::endl;
  }

  auto [time, solution] = solver.solution().back();

  CHECK(time == nT*step);
  for (std::ptrdiff_t col{0ll}; col < solution.cols(); ++col)
    for (std::ptrdiff_t row{0ll}; row < solution.rows(); ++row)
    {
      INFO("" << "row: " << row << ", col: " << col << ", solution: " << solution(row, col));
      CHECK_THAT(solution(row, col), WithinRel(val, tol));
    }
}
