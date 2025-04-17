
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::CoordinateTypes;

TEST_CASE("LogsTest")
{
    auto grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};

    auto permeability_stencils{Logs::Factory::generate_permeability_StepProperty(grid_stencils)};

    auto porosity_stencils{Logs::Factory::generate_porosity_StepProperty(grid_stencils)};

    auto grid{AxesGrid<Z>{GridDual{grid_stencils}}};
    auto permeability{Permeability{StepProperty{permeability_stencils}, grid}};
    auto porosity{Porosity{StepProperty{porosity_stencils}, grid}};
    
}
