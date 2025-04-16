#include <Injector/Grids/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;

// Tests Cylinder grid, (r; z)
TEST_CASE("GridTest", "GeneralCoordinate")
{
    Factory::create_cartesian_grid_2D(2);
    Factory::create_cylinder_grid_2D(5);
}
