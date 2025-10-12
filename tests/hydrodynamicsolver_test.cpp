#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::CoordinateTypes;

const auto z_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
const auto r_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};
// hydrodynamic logs
const auto is_permeable_stencils{
    Logs::RawDataFactory::generate_is_permeable(z_stencils)};
auto is_perforated_stencils{
    Logs::RawDataFactory::generate_is_permeable(z_stencils)};
const auto porosity_stencils{
    Logs::RawDataFactory::generate_porosity(z_stencils, is_permeable_stencils)};
const auto permeability_stencils{
    Logs::RawDataFactory::generate_permeability(z_stencils, is_permeable_stencils)};
const auto medium_compressibility_stencils{
    Logs::RawDataFactory::generate_medium_compressibility(z_stencils, is_permeable_stencils)};
const auto ext_pressure_stencils{
    Logs::RawDataFactory::generate_ext_pressure(z_stencils, is_permeable_stencils)};
// heat logs
const auto solid_density_stencils{
    Logs::RawDataFactory::generate_solid_density(z_stencils)};
const auto solid_specific_heatcapacity_stencils{
    Logs::RawDataFactory::generate_solid_specific_heatcapacity(z_stencils)};
const auto solid_heatconductivity_stencils{
    Logs::RawDataFactory::generate_conductivity(z_stencils)};

TEST_CASE("HydrodynamicsSolverTest")
{
    auto it = std::ranges::find_if(
        is_perforated_stencils,
        [](RealType v)
        { return v == 1.0; });
    (*it) = 0.0;

    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
    const auto &grid_z{grid2D->first_coord()};

    const Logs::Rocks::CoreSampleLogs core_data{
        is_permeable_stencils,
        is_perforated_stencils,
        porosity_stencils,
        permeability_stencils,
        grid_z};

    const Logs::Hydrodynamics::BaseHydrodynamics
        base_hydrodynamics{
            core_data,
            medium_compressibility_stencils,
            ext_pressure_stencils,
            grid_z};

    const Properties::Rocks::RocksProps
        collector_field{
            base_hydrodynamics, grid2D};

    const Logs::Rocks::HeatLogs
        heat_logs{
            solid_density_stencils,
            solid_specific_heatcapacity_stencils,
            solid_heatconductivity_stencils,
            porosity_stencils,
            Phases::FluidFactory::create_water(1.0, 1.0),
            grid_z};

    const Properties::Rocks::HeatProps
        heat_props{
            heat_logs, grid2D};
}
