
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Grids/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::Problem;
using namespace GPN::EqSolver::State;

// Tests Cartesian grid
TEST_CASE("PhaseProperties", "Water")
{
    auto z_grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
    auto z_grid{ZGrid{GridDual{z_grid_stencils}}};
    auto r_grid{RGrid{
        GridDual{
            GridDualStencils{
                Grids::Factory::
                    generate_dual_grid_stencils_uniform(0, 1, 11)}}}};
#pragma region GRID_2D
    // generate 1D grids in every direction --- points of property jumps
    const cptr<StructuredCylinderGrid2DAxisymmetric>
        grid2D{std::make_shared<StructuredCylinderGrid2DAxisymmetric>(
            z_grid, r_grid)};
#pragma endregion

    InitialCondition
        init_cond{
            State2D::FillWithConst(
                *grid2D, 1.0)};
}
