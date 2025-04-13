#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include <src/Grid/Defines.h>
#include <src/Solver/SplittingMethod/Solver.hpp>
#include <src/Grid/UniformGrid2D.h>

#include <src/Grid/UniformGridFactory.hpp>

#include <src/ConcreteProblems/SimpleProblemFactory.hpp>
#include <src/ConcreteProblems/FactoryWithSource.hpp>

using namespace EqSolver;
using namespace EqSolver::Grid;
using namespace EqSolver::Properties;
using namespace EqSolver::ConcreteProblem;
using namespace EqSolver::SplittingMethod;

TEST(Solver, splitting_method_homogeneous)
{
  const double tol = 1E-8;

  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.11, 0.11})};

  SimpleProblemFactory factory{grid};

  std::shared_ptr<Properties::Fields>
      properties{
          std::make_shared<Properties::Fields>(
              grid, factory)};

  Solver<Grid::UniformGrid2D> solver{
      grid, properties, factory};

  solver.advance(1);

  auto [time, state] = solver.solution().back();

  for (ptrdiff_t j = 0; j < (ptrdiff_t)grid.Y_nodes.size(); ++j)
    for (ptrdiff_t i = 0; i < (ptrdiff_t)grid.X_nodes.size(); ++i)
    {
      ASSERT_NEAR(state(i, j), 0.0, tol);
    }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}