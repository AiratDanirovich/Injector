#include <iostream>
#include <memory>

#include <catch2/catch_test_macros.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Grids/PhysicalField.hpp>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Solver/SplittingMethod/SplitX.hpp>

using namespace GPN;
using namespace GPN::EqSolver;
using namespace GPN::Grids;
using namespace GPN::EqSolver::SplittingMethod;

TEST_CASE("Solver", "splitX")
{
#pragma region GRID_2D
  // generate 1D grids in every direction --- points of property jumps
  StructuredCylinderGrid2DAxisymmetric
      grid2D{
          ZGrid{
              GridDual{
                  GridDualStencils{
                      Grids::Factory::
                          generate_dual_grid_stencils_uniform(0, 2, 3)}}},
          RGrid{
              GridDual{
                  GridDualStencils{
                      Grids::Factory::
                          generate_dual_grid_stencils_uniform(0, 2, 3)}}}};
#pragma endregion
#pragma region HEAT-CONDUCTIVITY
  // generate heat conductivity field
  Properties::HeatConductivity conductivity_field{
      Logs::HeatConductivity{
          Logs::StepPropertyGrid{
              Logs::StepProperty{
                  Logs::Factory::generate_conductivity_StepProperty(
                      grid2D.first_coord.dual_stencils)},
              grid2D.first_coord}},
      grid2D};
#pragma endregion

  const double tol = 1E-8;
  SplitX splitx{conductivity_field, grid2D};

  for (const auto &m : splitx.LaplaceTerms())
  {
    REQUIRE(m.rows() == m.cols());
    {
      auto col{0ll};
      CHECK(m.coeff(col, col) + m.coeff(col, col + 1ll) == 0.0);
    }
    for (auto col{1ll}; col < m.cols() - 1; ++col)
    {
      CHECK(m.coeff(col, col - 1ll) + m.coeff(col, col) + m.coeff(col, col + 1ll) == 0.0);
      CHECK(m.coeff(col, col) > tol);
    }
    {
      auto col{m.cols() - 1ll};
      CHECK(m.coeff(col, col - 1ll) + m.coeff(col, col) == 0.0);
    }
  }
}