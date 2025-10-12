#include <iostream>

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/transfer_to_vector.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace std;
using namespace GPN;
using namespace Catch;
using namespace Catch::Matchers;

using VR = std::vector<RealType>;

VR generate_stencils(const VR &steps, const RealType top = 0.0)
{
    VR out;
    out.reserve(steps.size() + 1ull);
    out.push_back(top);

    for (auto i{0ull}; i < steps.size(); ++i)
        out.push_back(out.back() + steps[i]);
    return out;
}

/*input data*/
// z-grid data
const auto nLayers{5ull};
const auto z_grid_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};

const RealType z_top{0.0};
const auto geotherma_nodes{transfer_to_vector(
    Logs::RawDataFactory::generate_geotherma_nodes(z_top - 1.0, z_top + 2.0, 31))};
const RealType ref_node{-2.0}, ref_val{0.0}, slope{1.5};
const auto geotherma_vals{
    Logs::RawDataFactory::generate_geotherma_vals_linear(
        geotherma_nodes, ref_node, ref_val, slope)};

TEST_CASE("GeothermaTest")
{
    {
        const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(z_grid_stencils)};
        const auto geotherma{Logs::GeothermaFactory::create(
            geotherma_nodes,
            geotherma_vals,
            z_top,
            grid)};

        const auto &mesh{grid.mesh_nodes};
        CHECK(geotherma.size() == grid.mesh_nodes.size());
        for (auto i{0ll}; i < geotherma.size(); ++i)
            CHECK_THAT(ref_val + slope * (mesh(i) - ref_node), WithinRel(geotherma(i), 1E-10));
    }
    {
        const RealType z_minor_step{3.0};
        LogValuesContainer is_permeable_stencils(9ll);
        is_permeable_stencils <<0.0, 0.0,  1.0,  0.0, 1.0,   0.0, 1.0, 1.0,  0.0;
        const VR thickness{7, 3, 3, 15, 3, 3, 3, 3, 35};
        const auto z_stencils{generate_stencils(thickness)};
        GPN::Grids::RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
        const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(refiner, z_stencils)};

        const VR geotherma_nodes{0.0, 60.0, 300, 2000, 2503, 2508, 3000};
        const VR geotherma_vals{293, 273.0, 283, 323, 348, 348, 364};
        const auto geotherma{Logs::GeothermaFactory::create(
            geotherma_nodes,
            geotherma_vals,
            z_top,
            grid)};

            Eigen::MatrixX<RealType> out(2, grid.mesh_nodes.size());
            out.row(0ll) = grid.mesh_nodes.transpose();
            out.row(1ll) =  geotherma.log_vals.transpose();

            cout << out << endl;
    }
}
