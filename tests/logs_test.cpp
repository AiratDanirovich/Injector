
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include <Injector/Model/Collector.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;

/*input data*/
// z-grid data
const auto nLayers{5ull};
const auto grid_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};
// hydrodynamic logs
const auto is_permeable_stencils{
    Logs::RawDataFactory::generate_is_permeable(grid_stencils)};
const auto porosity_stencils{
    Logs::RawDataFactory::generate_porosity(grid_stencils, is_permeable_stencils)};
const auto permeability_stencils{
    Logs::RawDataFactory::generate_permeability(grid_stencils, is_permeable_stencils)};

const auto ext_pressure_stencils{
    Logs::RawDataFactory::generate_ext_pressure(grid_stencils, is_permeable_stencils)};
const auto skin_stencils{
    Logs::RawDataFactory::generate_skin_StepProperty(grid_stencils, is_permeable_stencils)};
// heat logs
const auto solid_density_stencils{
    Logs::RawDataFactory::generate_solid_density(grid_stencils)};
const auto solid_specific_heatcapacity_stencils{
    Logs::RawDataFactory::generate_solid_specific_heatcapacity(grid_stencils)};
const auto heatconductivity_stencils{
    Logs::RawDataFactory::generate_conductivity_StepProperty(grid_stencils)};



TEST_CASE("LogsTest")
{
    const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(grid_stencils)};

    const Logs::Rocks::CoreSampleLogs core_data{
        is_permeable_stencils,
        porosity_stencils,
        permeability_stencils,
        grid};

    const Logs::Hydrodynamics::Hydrodynamics hydrodynamics{
        is_permeable_stencils,
        ext_pressure_stencils,
        skin_stencils, grid};
        
    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        heatconductivity_stencils,
        grid};
}
