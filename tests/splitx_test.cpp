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
  const double tol = 1E-8;
  #pragma region GRID_2D
  // generate 1D grids in every direction --- points of property jumps
  StructuredCylinderGrid2DAxisymmetric
      grid2D{
          ZGrid{
              GridDual{
                  GridDualStencils{
                      Grids::Factory::
                          generate_dual_grid_stencils_uniform(0, 1, 5)}}},
          RGrid{
              GridDual{
                  GridDualStencils{
                      Grids::Factory::
                          generate_dual_grid_stencils_uniform(0, 1, 11)}}}};
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

//  SplitX splitx{properties, grid};
}