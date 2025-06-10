#include <iostream>

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace GPN;
using namespace std;
using namespace Catch;
using namespace Catch::Matchers;

using VR = std::vector<RealType>;
LogValuesContainer transfer_to_eigen(const VR &data)
{
    LogValuesContainer out(data.size());
    for (auto i{0ull}; i < data.size(); ++i)
        out(i) = data[i];
    return out;
}

/*input data*/
// z-grid data
const auto nLayers{5ull};
const auto z_grid_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};

const RealType z_top{0.0};
const auto geotherma_nodes{
    Logs::RawDataFactory::generate_geotherma_nodes(z_top - 1.0, z_top + 2.0, 31)};
const RealType ref_node{-2.0}, ref_val{0.0}, slope{1.5};
const auto geotherma_vals{
    Logs::RawDataFactory::generate_geotherma_vals_linear(
        geotherma_nodes, ref_node, ref_val, slope)};

TEST_CASE("FieldsTest")
{
    const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(z_grid_stencils)};
    const auto geotherma{Logs::GeothermaFactory::create(
        geotherma_nodes,
        geotherma_vals,
        z_top,
        grid)};

        const auto& mesh{grid.mesh_nodes};
        CHECK(geotherma.size() == grid.mesh_nodes.size());
        for(auto i{0ll}; i < geotherma.size(); ++ i)
            CHECK_THAT(ref_val + slope*(mesh(i) - ref_node), WithinRel( geotherma(i), 1E-10));
}
