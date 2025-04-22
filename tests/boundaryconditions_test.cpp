
#include <memory>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Grids/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::BoundaryConditions;

struct BCFunctor : public GPN::BoundaryConditions::BCFunctorBase
{
    RealType operator()(
        RealType x, RealType y, RealType t) const override
    {
        return 1.0;
    }
};

// Tests Cartesian grid
TEST_CASE("PhaseProperties", "Water")
{
    auto z_grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
    auto z_grid{ZGrid{GridDual{z_grid_stencils}}};
#pragma region GRID_2D
    // generate 1D grids in every direction --- points of property jumps
    StructuredCylinderGrid2DAxisymmetric
        grid2D{z_grid,
               RGrid{
                   GridDual{
                       GridDualStencils{
                           Grids::Factory::
                               generate_dual_grid_stencils_uniform(0, 1, 11)}}}};
#pragma endregion

GPN::BoundaryConditions::BoundaryConditions bc{grid2D, std::make_shared<BCFunctor>(BCFunctor{})};
}
