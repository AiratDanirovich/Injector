#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;

// Tests Cylinder grid, (r; z)
TEST_CASE("GridTest", "GeneralCoordinate")
{
    const auto z_dual_size{3ll};
    const auto r_dual_size{5ll};
    auto z_stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, z_dual_size)};
    auto r_stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, r_dual_size)};
    {
        auto nodes{GridDual{z_stencils}};

        auto x_grid{
            AxesGrid<CoordinateTypes::X>{nodes}};

        auto y_grid{
            AxesGrid<CoordinateTypes::Y>{nodes}};

        StructuredXYGrid2D result{x_grid, y_grid};
    }

    {
        auto nodes{GridDual{r_stencils}};

        auto z_grid{
            AxesGrid<CoordinateTypes::Z>{nodes}};

        auto r_grid{
            AxesGrid<CoordinateTypes::R_CylCoord>{nodes}};

        StructuredCylinderGrid2DAxisymmetric result{z_grid, r_grid};
    }
    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
}
