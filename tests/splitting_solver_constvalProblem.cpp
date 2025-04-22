
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
  const double tol = 1E-8;

  RealType val{1.0};

  Solver solver{Factory::make_solver(val)};

  RealType step{0.0001};
  solver.advance(step);

  auto [time, solution] = solver.solution().back();

  CHECK(time == step);
  for (std::ptrdiff_t col{0ll}; col < solution.cols(); ++col)
    for (std::ptrdiff_t row{0ll}; row < solution.rows(); ++row)
    {
      INFO("" << "row: " << row << ", col: " << col << ", solution: " << solution(row, col));
      CHECK_THAT( solution(row, col), WithinRel(val, 1e-11) );
    }
}
