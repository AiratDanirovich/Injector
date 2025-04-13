#include <iostream>
#include <memory>

#include <catch2/catch_test_macros.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids.h>
#include <Injector/Grids/Factory/Factory.h>
#include <Injector/Solver/SplittingMethod/SplitX.hpp>

// #include <Injector/Grid/UniformGridFactory.hpp>

using namespace GPN;
using namespace GPN::EqSolver;
// using namespace GPN::Grids;
using namespace GPN::EqSolver::Properties;
using namespace GPN::EqSolver::SplittingMethod;

struct Capacity
{
  RealType operator()(RealType x, RealType y) const
  {
    return 1.0;
  }
};

struct Conductivity
{
  RealType operator()(RealType x, RealType y) const
  {
    return 1.0;
  }
};

struct Source
{
  RealType operator()(RealType x, RealType y) const
  {
    return 1.0;
  }
};

TEST_CASE("Solver", "splitX")
{
  const double tol = 1E-8;

  auto grid{
      Grids::Factory::create_grid_2D(3)
  };

  auto properties{
          std::make_shared<Properties::Fields>(
              grid,
              Capacity{},
              Conductivity{},
              Source{})
  };

  SplitX splitx{properties, grid};
}