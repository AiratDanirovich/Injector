#include <iostream>

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace std;

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

const auto geotherma_nodes{
    Logs::RawDataFactory::generate_geotherma_nodes(-1.0, 2.0, 31)};
const RealType ref_node{-2.0}, ref_val{0.0}, slope{1.5};
const auto geotherma_vals{
    Logs::RawDataFactory::generate_geotherma_vals_linear(
        geotherma_nodes, ref_node, ref_val, slope)};

TEST_CASE("FieldsTest")
{
    {
        const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(z_grid_stencils)};

        const auto &z_mesh{grid.mesh_nodes};
    }
    {
    }
}
