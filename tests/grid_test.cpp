#include <Injector/Grids/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;

// Tests Cylinder grid, (r; z)
TEST_CASE("GridTest", "GeneralCoordinate")
{
    {
        //    Factory::create_cartesian_grid_2D(2);
        auto stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, 2)};

        auto nodes{SpacialGrid{stencils}};

        auto x_grid{
            AxesGrid<CoordinateTypes::X>{nodes}};

        auto y_grid{
            AxesGrid<CoordinateTypes::Y>{nodes}};

        StructuredXYGrid2D result{x_grid, y_grid};
    }

    // Factory::create_cylinder_grid_2D(5);

    {
        auto stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};

        auto nodes{GridDual{stencils}};

        auto z_grid{
            AxesGrid<CoordinateTypes::Z>{nodes}};

        auto r_grid{
            AxesGrid<CoordinateTypes::R_CylCoord>{nodes}};

        StructuredCylinderGrid2D result{z_grid, r_grid};
    }
}
