
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>
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
                    Logs::Factory::generate_conductivity_StepProperty(
                        grid2D->first_coord.dual_stencils)},
                grid2D->first_coord}},
        grid2D};
#pragma endregion
#pragma region BASE-SPLIT
    BaseSplit base_split{
        conductivity_field.face_vals_axes2,
        grid2D->first_coord.mesh_size(),
        grid2D->second_coord.mesh_size()};
#pragma endregion
}
