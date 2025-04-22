
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Grids/PhysicalField.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <Injector/Solver/SplittingMethod/BaseSplit.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::EqSolver::SplittingMethod;

// Tests Cylinder grid, (r; z)
TEST_CASE("BaseSplitTest")
{
#pragma region GRID_2D
    // generate 1D grids in every direction --- points of property jumps
    const cptr<StructuredCylinderGrid2DAxisymmetric> grid2D{
        std::make_shared<StructuredCylinderGrid2DAxisymmetric>(
            Grids::Factory::create_cylinder_grid_2D(
                Box{Segment{0, 1}, Segment{0, 1}}, 5, 11))};
#pragma endregion
#pragma region HEAT-CONDUCTIVITY
    // generate heat conductivity field
    Properties::HeatConductivity conductivity_field{
        Logs::HeatConductivity{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    Logs::Factory::generate_conductivity_StepProperty(
                        grid2D->first_coord.dual_stencils)},
                grid2D->first_coord}},
        grid2D};
#pragma endregion
#pragma region BASE-SPLIT
    BaseSplit base_split{
        conductivity_field.face_vals_axes2,
        grid2D->first_coord.size(),
        grid2D->second_coord.size()};
#pragma endregion
}
