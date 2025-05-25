#include <memory>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>

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
    const auto grid2D{CylinderGridFactory::create(z_stencils, r_stencils)};
    BoundaryConditions::BoundaryConditions bc{*grid2D, std::make_shared<BCFunctor>(BCFunctor{})};
}
