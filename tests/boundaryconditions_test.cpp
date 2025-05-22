
#include <memory>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Grids/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;

struct BCFunctor : public BoundaryConditions::BCFunctorBase
{
    RealType operator()(
        RealType x, RealType y, RealType t) const override
    {
        return 1.0;
    }
};

auto z_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
auto r_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};

// Tests Cartesian grid
TEST_CASE("PhaseProperties", "Water")
{
    const CylinderGridFactory grid_factory{z_stencils, r_stencils};
    BoundaryConditions::BoundaryConditions bc{*grid_factory.grid(), std::make_shared<BCFunctor>(BCFunctor{})};
}
