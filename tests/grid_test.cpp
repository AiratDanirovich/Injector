#include <iostream>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Grids/GridRefiners.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::Grids;

using namespace std;

using VR = std::vector<RealType>;
LogValuesContainer transfer_to_eigen(const VR &data)
{
    LogValuesContainer out(data.size());
    for (auto i{0ull}; i < data.size(); ++i)
        out(i) = data[i];
    return out;
}

// Tests Cylinder grid, (r; z)
TEST_CASE("GridTest", "GeneralCoordinate")
{
    const RealType tol = 1e-11;

    const auto z_dual_size{3ll};
    const auto r_dual_size{15ll};
    auto z_stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, z_dual_size)};
    auto r_stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, r_dual_size)};
    {
        cout << "z_stencils:\n"
             << transfer_to_eigen(z_stencils).transpose() << endl;

        auto nodes{GridDual{z_stencils, CoordinateTypes::Z{}}};
        cout << "z dual nodes:\n"
             << nodes.dual_nodes.transpose() << endl;
        cout << "z mesh nodes:\n"
             << nodes.mesh_nodes.transpose() << endl;
#pragma region
        CHECK(nodes.dual_size() == z_dual_size);
        for (auto i{0ll}; i < (ptrdiff_t)z_stencils.size(); ++i)
            CHECK(nodes.dual_nodes(i) == z_stencils[i]);
#pragma endregion
#pragma region
        CHECK(nodes.mesh_size() == z_dual_size - 1ll);
        CHECK_THAT(nodes.mesh_nodes(0ll), WithinAbs(z_stencils[0], tol));
        for (auto i{1ll}; i < nodes.mesh_size() - 1ll; ++i)
        {
            CHECK(nodes.mesh_nodes(i) == (z_stencils[i] + z_stencils[i + 1ll]) / 2.0);
        }
        CHECK_THAT(nodes.mesh_nodes.tail(1ll)(0ll), WithinAbs(z_stencils.back(), tol));
#pragma endregion
        for (auto i{0ll}; i < nodes.dual_steps.size(); ++i)
            CHECK(nodes.dual_steps(i) == z_stencils[i + 1] - z_stencils[i]);

        auto x_grid{
            AxesGrid<CoordinateTypes::X>{nodes}};

        for (auto i{0ll}; i < x_grid.mesh_steps.size(); ++i)
            CHECK(x_grid.mesh_steps(i) == x_grid.mesh_nodes(i + 1ll) - x_grid.mesh_nodes(i));
        for (auto i{0ll}; i < x_grid.control_volumes.size(); ++i)
            CHECK(x_grid.control_volumes(i) == z_stencils[i + 1] - z_stencils[i]);

        auto y_grid{
            AxesGrid<CoordinateTypes::Y>{nodes}};

        StructuredXYGrid2D result{x_grid, y_grid};
        for (auto col{0ll}; col < y_grid.mesh_nodes.size(); ++col)
        {
            for (auto row{0ll}; row < x_grid.mesh_nodes.size(); ++row)
                CHECK(result.volume(row, col) == result.second_coord.control_volumes(col) * result.first_coord.control_volumes(row));
        }
    }

    {
        auto z_nodes{GridDual{z_stencils, CoordinateTypes::Z{}}};
        auto r_nodes{GridDual{r_stencils, CoordinateTypes::R_CylCoord{}}};

        cout << "r_stencils:\n"
             << transfer_to_eigen(r_stencils).transpose() << endl;

        auto z_grid{
            AxesGrid<CoordinateTypes::Z>{z_nodes}};

        auto r_grid{
            AxesGrid<CoordinateTypes::R_CylCoord>{r_nodes}};
        for (auto i{0ll}; i < r_grid.mesh_steps.size(); ++i)
            CHECK(r_grid.mesh_steps(i) == r_grid.mesh_nodes(i + 1ll) - r_grid.mesh_nodes(i));
        for (auto i{0ll}; i < r_grid.control_volumes.size(); ++i)
            CHECK(r_grid.control_volumes(i) == (r_stencils[i + 1] * r_stencils[i + 1] - r_stencils[i] * r_stencils[i]) / 2.0);

        StructuredCylinderGrid2DAxisymmetric result{z_grid, r_grid};
        for (auto col{0ll}; col < r_grid.mesh_nodes.size(); ++col)
        {
            for (auto row{0ll}; row < z_grid.mesh_nodes.size(); ++row)
                CHECK_THAT(result.volume(row, col),
                           WithinRel(
                               StructuredCylinderGrid2DAxisymmetric::TwoPI() *
                                   result.second_coord.control_volumes(col) *
                                   result.first_coord.control_volumes(row),
                               tol));
        }
    }
    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};


LogValuesContainer is_permeable(z_dual_size - 1ll);
is_permeable << 0.0 , 0.0;

RefinerVerticle refiner{0.01, is_permeable};

cout << refiner.refine(z_stencils);


}
