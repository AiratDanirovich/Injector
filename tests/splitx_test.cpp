#include <iostream>
#include <memory>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Solver/SplittingMethod/SplitX.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::EqSolver;
using namespace GPN::Grids;
using namespace GPN::EqSolver::SplittingMethod;

TEST_CASE("Solver", "splitX")
{
#pragma region GRID_2D
  // generate 1D grids in every direction --- points of property jumps
  const auto grid2D{
          Grids::Factory::create_cylinder_grid_2D_ptr(
              Box{Segment{0, 1}, Segment{0, 1}}, 5, 11)};
#pragma endregion
#pragma region HEAT-CONDUCTIVITY
  // generate heat conductivity field
  Properties::HeatConductivity conductivity_field{
      Logs::HeatConductivity{
          Logs::StepPropertyGrid{
              Logs::StepProperty{
                  Logs::RawDataFactory::generate_conductivity_StepProperty(
                      grid2D->first_coord.dual_stencils)},
              grid2D->first_coord}},
      grid2D};
#pragma endregion

  const double tol = 1E-15;
  SplitX splitx{conductivity_field, grid2D};

  for (const auto &m : splitx.LaplaceTerms())
  {
    REQUIRE(m.rows() == m.cols());
    {
      auto col{0ll};
      INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << 0.0 << ", r: " << m.coeff(col, col+1ll));
      CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col + 1ll), tol));
    }
    for (auto col{1ll}; col < m.cols() - 1; ++col)
    {
      INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " <<  m.coeff(col, col-1ll) << ", r: " << m.coeff(col, col+1ll));
      CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll) + m.coeff(col, col + 1ll), tol));
      INFO("" << "col: " << col << ", c: " << m.coeff(col, col));
      CHECK(m.coeff(col, col) > tol);
    }
    {
      auto col{m.cols() - 1ll};
      INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " <<  m.coeff(col, col-1ll) << ", r: " << 0.0);
      CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll), tol));
    }
  }
}