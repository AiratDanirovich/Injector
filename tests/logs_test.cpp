
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include <Injector/Model/Collector.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;

/*input data*/
// z-grid data
const auto nLayers{5ull};
const auto grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};
// hydrodynamic logs
const auto is_permeable_stencils{Logs::RawdataFactory::generate_is_permeable_StepProperty(grid_stencils)};
const auto porosity_stencils{Logs::RawdataFactory::generate_porosity_StepProperty(grid_stencils, is_permeable_stencils)};
const auto permeability_stencils{Logs::RawdataFactory::generate_permeability_StepProperty(grid_stencils, is_permeable_stencils)};

TEST_CASE("LogsTest")
{
    const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(grid_stencils)};
    const Logs::Rocks::CoreSampleData core_data{
        is_permeable_stencils,
        porosity_stencils,
        permeability_stencils,
        grid};
}
