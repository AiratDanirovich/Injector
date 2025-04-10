// #include <iostream>

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids.h>
#include <Injector/Grids/Factory/Factory.h>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
// using namespace GPN::Grids::Factory;

// Tests Cartesian grid
TEST_CASE("GridTest", "CartesianCoordinate")
{
    auto stencils{Factory::generate_z_stencils_uniform(0, 1 , 3)};

    auto nodes{GridNodes{stencils}};

    auto z_grid{Grid1D<CoordinateTypes::Z>{nodes}};

    auto r_grid{Grid1D<CoordinateTypes::RadialCylinderCoordinate>{nodes}};
}

// TEST_CASE("GridTest", "RadialCylinderCoordinate")
// {
    
// }

// int main(int argc, char **argv) 
// {
//     testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }