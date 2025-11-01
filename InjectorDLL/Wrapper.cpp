#include <vector>
#include <cmath>
#include <memory>
#include <numbers>
#include <fstream>
#include <iostream>
#include <filesystem>

#include <InjectorDLL/Wrapper.h>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Map/Grids2DMap.hpp>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/History/History.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluid.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleRatesFactory.hpp>
#include <Injector/Model/Heat/HeatBoundaryConditions.hpp>
#include <Injector/Model/Completion.hpp>
#include <Injector/Model/ExtrudedCasingFactory.hpp>

#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/FlowField.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include <tests/includes/get_completion.hpp>
#include <tests/includes/transfer_to_eigen.hpp>
#include <tests/includes/make_r_stencils.hpp>
#include <tests/includes/set_is_permeable_stencils.hpp>
#include <tests/includes/IC_BC.hpp>
#include <tests/includes/make_water.hpp>

#include <Eigen/Core>

#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

using namespace GPN;
using namespace GPN::CrossFlow;
using namespace GPN::Phases;
using namespace GPN::Completion;
using namespace GPN::Hydrodynamic;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::FullImplicit;

namespace fs = std::filesystem;

Wrapper::Wrapper(
    // fluid params in SI
    const RealType density,                 // kg/(m^3)
    const RealType capacity,                // J/(kg*K) /* specific heat capacity */
    const RealType viscosity,               // Pa*s
    const RealType heat_conductivity_fluid, // Watt/(m*K)
    const RealType joule_thomson,           // K/bar
    // grid
    const RealType rMin,         // m /* typically would be zero */
    const RealType rMax,         // m
    const RealType q,            // --, q >= 1.0 /* step increment factor */
    const RealType r_max_step,   // m /* maximum allowed step in radial direction */
    const RealType z_minor_step, // m, /*maximum step within impermeable layers*/
    // eight vectors of the same size
    // values are in SI
    const VR &thickness,                     // meter
    const VR &ext_pressure,                  // bar
    const VR &medium_compressibility,        // 1/Pa
    const VR &solid_heatconductivity,        // Watt/(m*K)
    const VR &porosity,                      // 0.0 < porosity <= 1.0, --
    const VR &permeability_stencils,         // m^2
    const VR &RFP_weights_stencils,          // -- /*rate distribution between layers*/
    const VR &is_perforated,                 // {0, 1}, --
    const VR &from_coords,                   // coordinates of column corrosion, m
    const std::vector<ptrdiff_t> &to_layers, // -- /* ids of layers accepting the cross flow */
    const VR &solid_density,                 // kg/(m^3)
    const VR &solid_specific_heatcapacity,   // J/(kg*K)
    // geotherma
    const RealType z_top,      // m, /* z-coordinate of the top */
    const VR &geotherma_nodes, // m, /* nodes for geotherma interpolation */
    const VR &geotherma_vals,  // K, /* reference vals for interpolation */
    // temporal grid
    const RealType start_time,   // start time in seconds
    const VR &time_intervals,    // intervals of const rates)
    const RealType t_minor_step, // time step used for numerical integration
    // well
    const VR &well_rates,         // ~1.1E-3 m^3/s
    const VR &inlet_temperatures, // K
    // casing
    // {fluid, tube, annulus, column, cement}
    const json &data)
{
    // adapt stl container to Eigen container
    //    LogValuesContainer is_permeable_stencils(is_permeable.size());
    //    std::copy(is_permeable.cbegin(), is_permeable.cend(), is_permeable_stencils.begin());
    LogValuesContainer is_perforated_stencils(is_perforated.size());
    std::copy(
        is_perforated.cbegin(),
        is_perforated.cend(),
        is_perforated_stencils.begin());
    LogValuesContainer ext_pressure_stencils{transfer_to_eigen(ext_pressure, 1e5)};
    LogValuesContainer medium_compressibility_stencils(medium_compressibility.size());
    std::copy(
        medium_compressibility.begin(),
        medium_compressibility.end(),
        medium_compressibility_stencils.begin());
    LogValuesContainer solid_density_stencils(solid_density.size());
    std::copy(
        solid_density.begin(),
        solid_density.end(),
        solid_density_stencils.begin());
    LogValuesContainer solid_specific_heatcapacity_stencils(
        solid_specific_heatcapacity.size());
    std::copy(
        solid_specific_heatcapacity.begin(),
        solid_specific_heatcapacity.end(),
        solid_specific_heatcapacity_stencils.begin());
    LogValuesContainer porosity_stencils(porosity.size());
    std::copy(
        porosity.begin(), porosity.end(),
        porosity_stencils.begin());
    LogValuesContainer solid_heatconductivity_stencils(
        solid_heatconductivity.size());
    std::copy(
        solid_heatconductivity.begin(),
        solid_heatconductivity.end(),
        solid_heatconductivity_stencils.begin());

    const auto is_permeable_stencils{
        set_is_permeable_stencils(
            is_perforated_stencils, to_layers)};

    // z-refiner
    const auto z_stencils{Grids::Factory::generate_dual_grid_stencils_from_steps(
        0.0, thickness)};
    GPN::Grids::RefinerVerticle z_refiner{z_minor_step, is_permeable_stencils};

    const auto temp_grid_z{
        Grids::Factory::create_axes<CoordinateTypes::Z>(
            z_refiner, z_stencils)};
    const Casing<VarRing> completion{get_completion(data, temp_grid_z)};

    //    const Casing completion{make_completion(casing_data)};

    // r_stencils
    const VR r_stencils{
        WellHoles{WellHolesFactory::create(completion)}.get_stencils(
            data["grid"]["r_start"],
            data["grid"]["r_end"])};
    // r-refiner
    const Grids::AbstractRefinerRadial *r_refiner{
        make_r_refiner(data)};

    const auto r_nodes{r_refiner->refine(r_stencils)};
    // the grid itself
    const auto grid2D{
        Grids::CylinderGridFactory::create(
            z_refiner, z_stencils,
            r_nodes)};
    const auto &grid_r{grid2D->second_coord()};
    const auto &grid_z{grid2D->first_coord()};

    const cptr<Grids::CylinderGridRock> grid2D_rocks{
        make_shared<Grids::CylinderGridRock>(grid2D)};
    constexpr auto left_margin{3ll};
    const auto &grid_rocks_z{grid2D_rocks->first_coord()};
    const auto &grid_rocks_r{grid2D_rocks->second_coord()};

    // make fluid
    const PhasePropertiesJT water{make_water(data)};

    const ExtrudedCasing extr_completion{
        VarExtrudedCasingFactory::create(completion)};
    // collector
    const Logs::Rocks::CoreSampleLogs core_logs{
        is_permeable_stencils,
        is_perforated_stencils,
        porosity_stencils,
        permeability_stencils,
        grid_z};

    const Logs::Hydrodynamics::BaseHydrodynamics
        base_hydrodynamics{
            Logs::Rocks::CoreSampleLogs{
                is_permeable_stencils,
                is_perforated_stencils,
                porosity_stencils,
                permeability_stencils,
                grid_z},
            medium_compressibility_stencils,
            ext_pressure_stencils};

    const Properties::Rocks::RocksProps
        rock_field_props{
            base_hydrodynamics,
            water,
            grid2D_rocks};

    const auto RFP_weights{
        Logs::RFPFactory::create_from_container(
            RFP_weights_stencils,
            core_logs.is_permeable)};
    const CrossFlows cross_flows{
        RFP_weights, from_coords, to_layers};
    const auto WFP_weights{
        create_WFP(
            core_logs.is_perforated,
            RFP_weights,
            cross_flows)};

    const Well_CrossFlow well{RFP_weights, WFP_weights, cross_flows};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        solid_heatconductivity_stencils,
        porosity_stencils,
        water,
        grid2D->first_coord()};

    Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};
    heat_props.apply_well(extr_completion, well);

    FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};
    heat_face_props.apply_well(extr_completion, well);

    // history
    const ptr<History> history{make_shared<History>(
        HistoryFactory::createFixedRate(
            time_intervals, well_rates, inlet_temperatures))};
    // external pressure log
    const auto external_pressure{
        Logs::ExtPressureFactory::create(
            ext_pressure_stencils,
            is_permeable_stencils,
            grid_z)};
    // fluid model for the pressure field
    using FluidField_t =
        decltype(CompressibleFluidField{
            start_time,
            water,
            rock_field_props,
            well,
            history,
            grid2D_rocks});

    auto ptr_pressure_field{
        make_shared<FluidField_t>(
            start_time,
            water,
            rock_field_props,
            well,
            history,
            grid2D_rocks)};

    // rates field factory
    using RatesFactory_t =
        decltype(FaceProperties::CompressibleRatesFactory{
            ptr_pressure_field,
            grid2D_rocks, well, history, water});
    auto ptr_rates_factory{make_shared<RatesFactory_t>(
        ptr_pressure_field,
        grid2D_rocks, well, history, water)};
    // initial condition
    std::unique_ptr<const Logs::Geotherma> geotherma{
        make_unique<Logs::Geotherma>(
            Logs::GeothermaFactory::create(
                geotherma_nodes,
                geotherma_vals,
                z_top,
                grid2D->first_coord()))};

    const auto initial_state{GPN::ICFactory(start_time, grid2D, *geotherma)};
    // boundary conditions
    const GPN::Heat::HeatBC bc{
        grid2D,
        std::make_shared<GPN::FunctorBC<
            Well_CrossFlow,
            History,
            FluidField_t,
            RatesFactory_t>>(
            history, ptr_rates_factory, *geotherma, grid2D),
        ptr_rates_factory};
    // solver
    using Solver_t = decltype(Solver{
        heat_face_props.medium_heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        ptr_rates_factory, initial_state,
        bc, start_time});

    auto solver_ptr{std::make_shared<Solver_t>(
        heat_face_props.medium_heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        ptr_rates_factory, initial_state,
        bc, start_time)};

    const auto &solver{*solver_ptr};

    SolverManager solver_manager{history, solver_ptr};

    solver_manager.run(t_minor_step);

    const auto &[times, states] = solver.solution();
    this->time = times;
    
    const auto& [p_times, p_states] = ptr_rates_factory->solution;

    if (!fs::is_directory("output") || !fs::exists("output")) // Check if src folder exists
    {
        fs::create_directory("output"); // create src folder
    }
    if (!fs::is_directory("output/pressure") || !fs::exists("output/pressure")) // Check if src folder exists
    {
        fs::create_directory("output/pressure"); // create src folder
    }

    try
    {
        ifstream f("separators.json");
        json data = json::parse(f);

        const std::string sep = data["coeff_sep"];

        auto grid{grid2D->second_coord().mesh_nodes};
        grid(0ll) = grid_r.dual_nodes(1ll);

        const Eigen::IOFormat commaFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, sep, sep, "", "", "", "");
        for (auto z{0ll}, layer_id{0ll}; z < grid2D->first_coord().mesh_nodes.size(); ++z)
        {
            if (core_logs.is_permeable(z) == 1.0)
            {
                {
                    ofstream f{std::string{"output/layer_"} + std::to_string(layer_id) + std::string{".csv"}};

                    f << sep << sep
                      << grid.transpose().format(commaFmt) << '\n';
                    for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                    {
                        f << t << sep << times[t] << sep
                          << states[t].cur_state.row(z).format(commaFmt) << '\n';
                    }

                    f.close();
                }
                {
                    ofstream f{std::string{"output/pressure/layer_"} + std::to_string(layer_id) + std::string{".csv"}};

                    f << sep << sep
                      << grid.transpose().format(commaFmt) << '\n';
                    for (auto t{0ll}; t < (ptrdiff_t)p_times.size(); ++t)
                    {
                        f << t << sep << p_times[t] << sep
                          << p_states[t].cur_state.row(z).format(commaFmt) << '\n';
                    }

                    f.close();

                }
                ++layer_id;
            }
        }

        {
            ofstream f{std::string{"output/well_temperature.csv"}};
            f << sep << sep << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                f << t << sep << times[t] << sep << states[t].cur_state.col(0ll).format(commaFmt) << '\n';
            }
            f.close();
        }

        {
            ofstream f{std::string{"output/cement_temperature.csv"}};
            f << sep << sep << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                f << t << sep << times[t] << sep << states[t].cur_state.col(1ll).format(commaFmt) << '\n';
            }
            f.close();
        }

        // print grids
        {
            ofstream f{std::string{"output/z_grid.csv"}};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f.close();
        }
        {
            ofstream f{std::string{"output/r_grid.csv"}};
            f << grid_r.mesh_nodes.transpose().format(commaFmt) << '\n';
            f.close();
        }
        {
            ofstream f{std::string{"output/geotherma.csv"}};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << initial_state.cur_state.col(0ll).transpose().format(commaFmt) << '\n';
            f.close();
        }
        {
            ofstream f{std::string{"output/casing_radius_pos.csv"}};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << extr_completion.radial_node_position().transpose().format(commaFmt) << '\n';
            f.close();
        }
        {
            ofstream f{std::string{"output/flow-props.csv"}};
            const auto ring{MaterialType::Flow};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].density().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].specific_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].volumetric_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].heat_conductivity().transpose().format(commaFmt) << '\n';

            f << extr_completion[ring].inner_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].outer_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].thickness.transpose().format(commaFmt) << '\n';

            f.close();
        }
        {
            ofstream f{std::string{"output/tube-props.csv"}};
            const auto ring{MaterialType::Tube};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].density().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].specific_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].volumetric_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].heat_conductivity().transpose().format(commaFmt) << '\n';

            f << extr_completion[ring].inner_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].outer_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].thickness.transpose().format(commaFmt) << '\n';

            f.close();
        }
        {
            ofstream f{std::string{"output/annulus-props.csv"}};
            const auto ring{MaterialType::Annulus};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].density().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].specific_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].volumetric_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].heat_conductivity().transpose().format(commaFmt) << '\n';

            f << extr_completion[ring].inner_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].outer_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].thickness.transpose().format(commaFmt) << '\n';

            f.close();
        }
        {
            ofstream f{std::string{"output/column-props.csv"}};
            const auto ring{MaterialType::Column};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].density().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].specific_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].volumetric_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].heat_conductivity().transpose().format(commaFmt) << '\n';

            f << extr_completion[ring].inner_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].outer_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].thickness.transpose().format(commaFmt) << '\n';

            f.close();
        }
        {
            ofstream f{std::string{"output/cement-props.csv"}};
            const auto ring{MaterialType::Cement};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].density().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].specific_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].volumetric_heat_capacity().transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].heat_conductivity().transpose().format(commaFmt) << '\n';

            f << extr_completion[ring].inner_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].outer_radius.transpose().format(commaFmt) << '\n';
            f << extr_completion[ring].thickness.transpose().format(commaFmt) << '\n';

            f.close();
        }
        // {
        //     ofstream f{std::string{"output/data.txt"}};
        //     f << "top collector height: " << grid_z.mesh_nodes(well.top_collector_cell_id()) << " m" << endl;
        //     f.close();
        // }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
}

std::vector<RealType> Wrapper::get_times() const
{
    return time;
}
// std::vector<std::vector<RealType>> Wrapper::get_temps() const
// {
//     return t_radial_distribution;
// }

Wrapper::~Wrapper()
{
}