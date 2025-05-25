
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Grids/GridsFactory.hpp>

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
    auto r_grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};
    const auto grid2D{Grids::CylinderGridFactory::create(z_grid_stencils, r_grid_stencils)};

    InitialCondition
        init_cond{
            State2D::FillWithConst(
                *grid2D, 1.0)};
}
