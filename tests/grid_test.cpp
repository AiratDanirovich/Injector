#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids.h>
#include <Injector/Grids/Factory/Factory.h>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;

// Tests Cartesian grid
TEST_CASE("GridTest", "CartesianCoordinate")
{
    auto stencils{Factory::generate_stencils_uniform(0, 1 , 3)};

    auto nodes{GridNodes{stencils}};

    auto z_grid{
        Grid1D<CoordinateTypes::Z>{nodes}
    };

    auto r_grid{
        Grid1D<CoordinateTypes::RadialCylinderCoordinate>{nodes}
    };

    auto grid_2D{
        StructuredCylinderGrid2D{
            Grid1D<CoordinateTypes::Z>{nodes}, 
            Grid1D<CoordinateTypes::RadialCylinderCoordinate>{nodes}}
    };
}
