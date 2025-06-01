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
        CHECK(nodes.dual_size() == z_dual_size);
        for (auto i{0ll}; i < (ptrdiff_t)z_stencils.size(); ++i)
            CHECK(nodes.dual_nodes(i) == z_stencils[i]);

        CHECK(nodes.mesh_size() == z_dual_size - 1ll);
        CHECK(nodes.mesh_nodes(0ll) == z_stencils[0]);
        for (auto i{1ll}; i < (ptrdiff_t)z_stencils.size() - 1ll; ++i)
            CHECK(nodes.mesh_nodes(i) == (z_stencils[i] + z_stencils[i + 1ll]) / 2.0);
        CHECK(nodes.mesh_nodes.tail(1ll)(0ll) == z_stencils.back());

        for (auto i{0ll}; i < nodes.dual_steps.size(); ++i)
            CHECK(nodes.dual_steps(i) == z_stencils[i + 1] - z_stencils[i]);

        auto x_grid{
            AxesGrid<CoordinateTypes::X>{nodes}};

        for (auto i{0ll}; i < x_grid.mesh_steps.size(); ++i)
            CHECK(x_grid.mesh_steps(i) == z_stencils[i + 1] - z_stencils[i]);
        for (auto i{0ll}; i < x_grid.control_volumes.size(); ++i)
            CHECK(x_grid.control_volumes(i) == z_stencils[i + 1] - z_stencils[i]);

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
