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

TEST(Solver, splitting_method_uniform)
{
  const double tol = 1E-6;

  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.011, 0.011})};

  FactoryWithSource factory{grid};

  std::shared_ptr<Properties::Fields>
      properties{
          std::make_shared<Properties::Fields>(
              grid, factory)};

  Solver<Grid::UniformGrid2D> solver{
      grid, properties, factory};

  {
    // first step
    double tau = 0.0001;
    const double final_time = 0.01;
    double t = tau;
    for(; t < final_time + tol; t+=tau )
      solver.advance(tau);

    auto [time, state] = solver.solution().back();

double mmax = 0.0;

    for (ptrdiff_t j = 0; j < (ptrdiff_t)grid.Y_nodes.size(); ++j)
      for (ptrdiff_t i = 0; i < (ptrdiff_t)grid.X_nodes.size(); ++i)
      {
        // ASSERT_NEAR(state(i, j), 1.0, tol);
       //   EXPECT_NEAR(state(i, j), tau, tol) << " i = " << i << "; j = " << j << std::endl;

          mmax = (mmax < abs(state(i, j) - final_time)) ? abs(state(i, j) - final_time) : mmax;
      }

      EXPECT_NEAR(mmax, 0.0, final_time*5e-3) << "max_diff = " << mmax;
  }

  // {
  //   // second step
  //   solver.advance(1);

  //   auto [time, state] = solver.solution().back();

  //   for (ptrdiff_t j = 0; j < (ptrdiff_t)grid.Y_nodes.size(); ++j)
  //     for (ptrdiff_t i = 0; i < (ptrdiff_t)grid.X_nodes.size(); ++i)
  //     {
  //       EXPECT_NEAR(state(i, j), 2.0, tol);
  //     }
  // }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}