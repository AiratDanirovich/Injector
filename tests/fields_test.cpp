#include <iostream>

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

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
const auto grid_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};
// hydrodynamic logs
const auto is_permeable_stencils{
    Logs::RawDataFactory::generate_is_permeable(grid_stencils)};
const auto is_perforated_stencils{
    Logs::RawDataFactory::generate_is_permeable(grid_stencils)};
const auto porosity_stencils{
    Logs::RawDataFactory::generate_porosity(grid_stencils, is_permeable_stencils)};
const auto permeability_stencils{
    Logs::RawDataFactory::generate_permeability(grid_stencils, is_permeable_stencils)};

const auto ext_pressure_stencils{
    Logs::RawDataFactory::generate_ext_pressure(grid_stencils, is_permeable_stencils)};
const auto skin_stencils{
    Logs::RawDataFactory::generate_skin(grid_stencils, is_permeable_stencils)};
// heat logs
const auto solid_density_stencils{
    Logs::RawDataFactory::generate_solid_density(grid_stencils)};
const auto solid_specific_heatcapacity_stencils{
    Logs::RawDataFactory::generate_solid_specific_heatcapacity(grid_stencils)};
const auto heatconductivity_stencils{
    Logs::RawDataFactory::generate_conductivity(grid_stencils)};

TEST_CASE("FieldsTest")
{
    {
        const auto grid2D{
            Grids::CylinderGridFactory::create(grid_stencils, grid_stencils)};

        const auto &grid{grid2D->first_coord};

        const Logs::Rocks::CoreSampleLogs core_data{
            is_permeable_stencils,
            is_perforated_stencils,
            porosity_stencils,
            permeability_stencils,
            grid};

        const Logs::Hydrodynamics::Hydrodynamics hydrodynamics_logs{
            is_permeable_stencils,
            ext_pressure_stencils,
            skin_stencils, grid};

        const Logs::Rocks::HeatLogs heat_logs{
            solid_density_stencils,
            solid_specific_heatcapacity_stencils,
            heatconductivity_stencils,
            porosity_stencils,
            Phases::FluidFactory::create_water(1.0, 1.0),
            grid};

        const Properties::Rocks::Rocks collector_field{
            core_data, grid2D};

        const Properties::Rocks::HeatProps heat_props{
            heat_logs, grid2D};
    }
    {
        Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
        Grids::RefinerVerticle refiner{0.1, transfer_to_eigen(is_permeable_stencils)};

        auto temp = Grids::Factory::create_axes<CoordinateTypes::Z>(refiner, grid_stencils);

        cout << "is permeable:   " << transfer_to_eigen(is_permeable_stencils).transpose().format(CommaInitFmt) << endl;
        cout << "refined nodes:  " << temp.dual_nodes.transpose().format(CommaInitFmt) << endl;
        cout << "dual nodes:     " << temp.dual_stencils.dual_nodes.transpose().format(CommaInitFmt) << endl;

        const auto grid2D{
            Grids::CylinderGridFactory::create(refiner, grid_stencils, grid_stencils)};

        const auto &grid{grid2D->first_coord};

        const Logs::Rocks::CoreSampleLogs core_data{
            is_permeable_stencils,
            is_perforated_stencils,
            porosity_stencils,
            permeability_stencils,
            grid};

        cout << "refined is_permeable: " << core_data.is_permeable.log_vals.transpose().format(CommaInitFmt) << endl;
        cout << "refined permeability: " << core_data.permeability.log_vals.transpose().format(CommaInitFmt) << endl;
        cout << "refined porosity:     " << core_data.porosity.log_vals.transpose().format(CommaInitFmt) << endl;

        const Logs::Hydrodynamics::Hydrodynamics hydrodynamics_logs{
            is_permeable_stencils,
            ext_pressure_stencils,
            skin_stencils, grid};
            
        cout << "refined ext pressure: " << hydrodynamics_logs.ext_pressure.log_vals.transpose().format(CommaInitFmt) << endl;
        cout << "refined skin:         " << hydrodynamics_logs.skin.log_vals.transpose().format(CommaInitFmt) << endl;

        const Logs::Rocks::HeatLogs heat_logs{
            solid_density_stencils,
            solid_specific_heatcapacity_stencils,
            heatconductivity_stencils,
            porosity_stencils,
            Phases::FluidFactory::create_water(1.0, 1.0),
            grid};
            
        cout << "refined density:             " << heat_logs.solid_density.log_vals.transpose().format(CommaInitFmt) << endl;
        cout << "refined spec heat cap:       " << heat_logs.solid_specific_heatcapacity.log_vals.transpose().format(CommaInitFmt) << endl;
        cout << "refined medium vol heat cap: " << heat_logs.medium_vol_heatcapacity.log_vals.transpose().format(CommaInitFmt) << endl;
    }
}
