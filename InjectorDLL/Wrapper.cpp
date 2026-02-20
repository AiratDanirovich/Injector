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
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Well/WellReservoirFlowProfileControl.hpp>
#include <Injector/Model/Well/WellBottomHolePressureControl.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
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
#include <tests/includes/make_history.hpp>

#include <Eigen/Core>

#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

using namespace GPN;
using namespace GPN::CrossFlow;
using namespace GPN::Phases;
using namespace GPN::Completion;
using namespace GPN::Hydrodynamic;
//using namespace GPN::Wells::ResFlowProfileControl;
//using namespace GPN::Wells::BotHolePresControl;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::FullImplicit;

namespace fs = std::filesystem;

template <typename Well_t>
struct WrapperFactory
{
    static void choose_well(const json &data)
    {
        /*START*/
        // input parameters
        /*fluid*/
        RealType
            viscosity{data["fluid"]["viscosity"].get<RealType>()},
            density{data["fluid"]["density"].get<RealType>()},
            capacity{data["fluid"]["specific_heat_capacity"].get<RealType>()},
            heat_conductivity{data["fluid"]["heat_conductivity"].get<RealType>()},
            joule_thomson{data["fluid"]["joule_thomson"].get<RealType>()};
        /*collector*/
        const auto thickness{data["collector"]["thickness"].get<VR>()};
        // hydrodynamic logs
        const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
        const auto ext_pressure_stencils{transfer_to_eigen(data["collector"]["external_pressure"].get<VR>(), 1e5)};
        const auto medium_compressibility_stencils{transfer_to_eigen(data["collector"]["medium_compressibility"].get<VR>())};
        const auto solid_density_stencils{transfer_to_eigen(data["collector"]["solidDensity"].get<VR>())};
        const auto solid_specific_heatcapacity_stencils{transfer_to_eigen(data["collector"]["solidSpecificHeatCapacity"].get<VR>())};
        const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"].get<VR>())};
        const auto solid_heatconductivity_stencils{transfer_to_eigen(data["collector"]["heatConductivity"].get<VR>())};
        const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>())};
        const auto RFP_weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};
        const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
        const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};
        const auto is_permeable_stencils{
            set_is_permeable_stencils(
                is_perforated_stencils, to_layers)};
        // heat logs
        /*grid*/
        const RealType
            rMin{data["grid"]["r_start"].get<RealType>()},
            rMax{data["grid"]["r_end"].get<RealType>()},
            q{data["grid"]["r_log_grid"]["q"].get<RealType>()},
            r_max_step{data["grid"]["r_log_grid"]["r_max_step"].get<RealType>()},
            z_minor_step{data["grid"]["z_minor_step"].get<RealType>()}; // m
        /*history*/
        const RealType
            start_time{read_start_time(data)};
        const RealType t_minor_step{read_minor_step(data)};
        /*temperatures*/
        // const VR well_rates = data["history"]["dynamic"]["well_rate"].get<VR>(); // m^3/s
        // const VR inlet_temperatures = data["history"]["dynamic"]["inlet_temperature"].get<VR>();
        /*well*/
        //    const auto casing{parse_completion(data)};
        /*END*/

        const auto &data2 = data["collector"]["geotherma"]["interpolate"];
        const VR geotherma_nodes = data2["z_nodes"].get<VR>();
        const VR geotherma_vals = data2["t_vals"].get<VR>();
        const RealType z_top = data2["z_top"].get<RealType>();

        cout << "Simulation is started." << endl;
        cout << "Please wait..." << endl;
#pragma region
        // z-refiner
        const auto z_stencils{Grids::Factory::generate_dual_grid_stencils_from_steps(
            0.0, thickness)};
        GPN::Grids::RefinerVerticle z_refiner{z_minor_step, is_permeable_stencils};

        const auto temp_grid_z{
            Grids::Factory::create_axes<CoordinateTypes::Z>(
                z_refiner, z_stencils)};
        const auto casing{get_completion(data, temp_grid_z)};
        const Casing<VarRing> completion{casing};

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

        const Logs::RFP_weights RFP_w{
            Logs::RFPFactory::create_from_container<Logs::RFP_weights>(
                Logs::StepProperty{RFP_weights_stencils},
                core_logs.is_permeable)};
        const CrossFlows cross_flows{
            from_coords, to_layers, RFP_w, core_logs.is_perforated};

        // history
        const ptr<History> history{make_shared<History>(make_history(data))};

        const ptr<Well_t> well{
            std::make_shared<Well_t>(
                rock_field_props,
                cross_flows,
                history,
                water,
                grid2D_rocks)};

        const Logs::Rocks::HeatLogs heat_logs{
            solid_density_stencils,
            solid_specific_heatcapacity_stencils,
            solid_heatconductivity_stencils,
            porosity_stencils,
            water,
            grid2D->first_coord()};

        Properties::Rocks::HeatProps heat_props{
            heat_logs, grid2D};
        heat_props.apply_well(extr_completion);

        FaceProperties::Rocks::HeatFaceProps heat_face_props{
            heat_props, grid2D};
        heat_face_props.apply_well(extr_completion);

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
                *well,
                history,
                grid2D_rocks});

        auto ptr_pressure_field{
            make_shared<FluidField_t>(
                start_time,
                water,
                rock_field_props,
                *well,
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

        // saved temperature in computation domain
        const auto &[times, states] = solver.solution();

        // saved pressure in layers
        const auto &[p_times, p_states] = ptr_rates_factory->solution;
        // saved rates in layers
        const auto &[q_times, q_states] = ptr_rates_factory->rates;
        // saved data in well
        const auto &well_data{well->solution};
        const auto &w_times{well_data.times};

        if (!fs::is_directory("output") || !fs::exists("output")) // Check if src folder exists
        {
            fs::create_directory("output"); // create src folder
        }
        if (!fs::is_directory("output/pressure") || !fs::exists("output/pressure")) // Check if src folder exists
        {
            fs::create_directory("output/pressure"); // create src folder
        }
        if (!fs::is_directory("output/rate") || !fs::exists("output/rate")) // Check if src folder exists
        {
            fs::create_directory("output/rate"); // create src folder
        }
        if (!fs::is_directory("output/well") || !fs::exists("output/well")) // Check if src folder exists
        {
            fs::create_directory("output/well"); // create src folder
        }

        try
        {
            ifstream f("separators.json");
            json data = json::parse(f);

            const auto sep{data["coeff_sep"].get<std::string>()};

            // insert sandface in a proper position
            const auto &grid_r{grid2D->second_coord()};
            const auto size{grid_r.dual_size()};

            const Eigen::IOFormat commaFmt(Eigen::FullPrecision, Eigen::DontAlignCols, sep, "\n", "", "", "", "");
            for (auto z{0ll}, layer_id{0ll}; z < grid2D->first_coord().mesh_size(); ++z)
            {
                MeshNodesContainerT grid{MeshNodesContainerT::Zero(grid_r.dual_size())};
                grid.tail(grid2D->second_coord().mesh_size()) = grid_r.mesh_nodes;
                grid.head(3ll) = grid_r.mesh_nodes.head(3ll);
                grid(3ll) = completion.sandface_radius(z);

                if (core_logs.is_permeable(z) == 1.0)
                {
                    { // temperature in layers
                        ofstream f{std::string{"output/t_layer_"} + std::to_string(layer_id) + std::string{".csv"}};

                        // print time ids
                        f << sep;
                        for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                            f << sep << t;
                        f << '\n';

                        // print time moments
                        f << sep << sep << transfer_to_eigen(times).transpose().format(commaFmt) << '\n';

                        // prepare data to print
                        GridNodeValues2D out_tt{GridNodeValues2D::Zero(size, (ptrdiff_t)times.size() + 2ll)};
                        out_tt.col(0ll) = grid_r.dual_nodes.transpose();
                        out_tt.col(1ll) = grid.transpose();
                        for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                        {
                            MeshNodesContainer out_t{MeshNodesContainer::Zero(size)};
                            const auto T{states[t].cur_state.row(z).transpose()};
                            out_t.head(3ll) = T.head(3ll);
                            out_t.tail(size - 3ll) = T.tail(size - 3ll);
                            out_t(3ll) = (out_t(2ll) + out_t(4ll)) / 2.0;

                            out_tt.col(t + 2ll) = out_t;
                        }
                        // print prepared data
                        f << out_tt.format(commaFmt) << '\n';

                        f.close();
                    }
                    { // pressure in layers
                        ofstream f{std::string{"output/pressure/p_layer_"} + std::to_string(layer_id) + std::string{".csv"}};

                        // print time ids
                        f << sep;
                        for (auto t{0ll}; t < (ptrdiff_t)p_times.size(); ++t)
                            f << sep << t;
                        f << '\n';

                        // print time moments
                        f << sep << sep << transfer_to_eigen(p_times).transpose().format(commaFmt) << '\n';

                        // prepare data to print
                        GridNodeValues2D out_pp{GridNodeValues2D::Zero(size, (ptrdiff_t)times.size() + 2ll)};
                        out_pp.col(0ll) = grid_r.dual_nodes.transpose();
                        out_pp.col(1ll) = grid.transpose();
                        for (auto t{0ll}; t < (ptrdiff_t)p_times.size(); ++t)
                        {
                            MeshNodesContainer out_p{MeshNodesContainer::Zero(size)};
                            const auto P{p_states[t].cur_state.row(z).transpose()};
                            out_p.head(3ll) = P.head(3ll);
                            out_p.tail(size - 3ll) = P.tail(size - 3ll);
                            out_p(3ll) = out_p(2ll);

                            out_pp.col(t + 2ll) = out_p;
                        }
                        // print prepared data
                        f << out_pp.format(commaFmt) << '\n';
                        f.close();
                    }
                    
                    {   // rates in layers
                        ofstream f{std::string{"output/rate/q_layer_"} + std::to_string(layer_id) + std::string{".csv"}};

                        // print time ids
                        for (auto t{0ll}; t < (ptrdiff_t)q_times.size(); ++t)
                            f << sep << t;
                        f << '\n';

                        // print time moments
                        f << sep << transfer_to_eigen(q_times).transpose().format(commaFmt) << '\n';

                        // prepare data to print
                        const auto& grid_r{grid_rocks_r};
                        const auto size_r{grid_r.dual_size()};
                        GridNodeValues2D out_qq{GridNodeValues2D::Zero(size_r, (ptrdiff_t)q_times.size() + 1ll)};
                        out_qq.col(0ll) = grid_r.dual_nodes.transpose();
                        for (auto t{0ll}; t < (ptrdiff_t)q_times.size(); ++t)
                        {
                            const auto Q{q_states[t].cur_state.row(z).transpose()};
                            out_qq.col(t + 1ll) = Q;
                        }
                        // print prepared data
                        f << out_qq.format(commaFmt) << '\n';
                        f.close();
                    }
                    ++layer_id;
                }
            }

            // well rfp
            {
                const auto size_z{grid_rocks_z.mesh_size()};
                ofstream f{std::string{"output/well/rfp.csv"}};
                // print time ids
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    f << sep << t;
                f << '\n';

                // print time moments
                f << sep << transfer_to_eigen(w_times).transpose().format(commaFmt) << '\n';

                // prepare data to print
                GridNodeValues2D out_tt{GridNodeValues2D::Zero(size_z, (ptrdiff_t)w_times.size() + 1ll)};
                out_tt.col(0ll) = grid_rocks_z.mesh_nodes;
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    out_tt.col(t + 1ll) = well_data.rfp[t].cur_state;
                // print prepared data
                f << out_tt.format(commaFmt) << '\n';
                f.close();
            }
            // well wfp
            {
                const auto size_z{grid_rocks_z.mesh_size()};
                ofstream f{std::string{"output/well/wfp.csv"}};
                // print time ids
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    f << sep << t;
                f << '\n';

                // print time moments
                f << sep << transfer_to_eigen(w_times).transpose().format(commaFmt) << '\n';

                // prepare data to print
                GridNodeValues2D out_tt{GridNodeValues2D::Zero(size_z, (ptrdiff_t)w_times.size() + 1ll)};
                out_tt.col(0ll) = grid_rocks_z.mesh_nodes;
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    out_tt.col(t + 1ll) = well_data.wfp[t].cur_state;
                // print prepared data
                f << out_tt.format(commaFmt) << '\n';
                f.close();
            }
            // well verticle_cement_flow
            {
                const auto size_z{grid_rocks_z.dual_size()};
                ofstream f{std::string{"output/well/verticle_cement_flow.csv"}};
                // print time ids
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    f << sep << t;
                f << '\n';

                // print time moments
                f << sep << transfer_to_eigen(w_times).transpose().format(commaFmt) << '\n';

                // prepare data to print
                GridNodeValues2D out_tt{GridNodeValues2D::Zero(size_z, (ptrdiff_t)w_times.size() + 1ll)};
                out_tt.col(0ll) = grid_rocks_z.dual_nodes;
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    out_tt.col(t + 1ll) = well_data.verticle_cement_flow[t].cur_state;
                // print prepared data
                f << out_tt.format(commaFmt) << '\n';
                f.close();
            }
            // well verticle_well_flow
            {
                const auto size_z{grid_rocks_z.dual_size()};
                ofstream f{std::string{"output/well/verticle_well_flow.csv"}};
                // print time ids
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    f << sep << t;
                f << '\n';

                // print time moments
                f << sep << transfer_to_eigen(w_times).transpose().format(commaFmt) << '\n';

                // prepare data to print
                GridNodeValues2D out_tt{GridNodeValues2D::Zero(size_z, (ptrdiff_t)w_times.size() + 1ll)};
                out_tt.col(0ll) = grid_rocks_z.dual_nodes;
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    out_tt.col(t + 1ll) = well_data.verticle_well_flow[t].cur_state;
                // print prepared data
                f << out_tt.format(commaFmt) << '\n';
                f.close();
            }
            
            // well_pressure
            {
                const auto size_z{grid_rocks_z.mesh_size()};
                ofstream f{std::string{"output/well/well_pressure.csv"}};
                // print time ids
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    f << sep << t;
                f << '\n';

                // print time moments
                f << sep << transfer_to_eigen(w_times).transpose().format(commaFmt) << '\n';

                // prepare data to print
                GridNodeValues2D out_tt{GridNodeValues2D::Zero(size_z, (ptrdiff_t)w_times.size() + 1ll)};
                out_tt.col(0ll) = grid_rocks_z.mesh_nodes;
                for (auto t{0ll}; t < (ptrdiff_t)w_times.size(); ++t)
                    out_tt.col(t + 1ll) = well_data.well_pressure[t].cur_state;
                // print prepared data
                f << out_tt.format(commaFmt) << '\n';
                f.close();
            }

            const auto size_z{grid_z.mesh_size()};
            {
                ofstream f{std::string{"output/well_temperature.csv"}};
                // print time ids
                f << sep;
                for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                    f << sep << t;
                f << '\n';

                // print time moments
                f << sep << sep << transfer_to_eigen(times).transpose().format(commaFmt) << '\n';

                // prepare data to print
                GridNodeValues2D out_tt{GridNodeValues2D::Zero(size_z, (ptrdiff_t)times.size() + 2ll)};
                out_tt.col(0ll) = grid_z.dual_nodes.tail(size_z);
                out_tt.col(1ll) = grid_z.mesh_nodes;
                for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                    out_tt.col(t + 2ll) = states[t].cur_state.col(0ll);
                // print prepared data
                f << out_tt.format(commaFmt) << '\n';
                f.close();
            }

            {
                ofstream f{std::string{"output/cement_temperature.csv"}};
                // print time ids
                f << sep;
                for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                    f << sep << t;
                f << '\n';
                // print time moments
                f << sep << sep << transfer_to_eigen(times).transpose().format(commaFmt) << '\n';
                // prepare data to print
                GridNodeValues2D out_tt{GridNodeValues2D::Zero(size_z, (ptrdiff_t)times.size() + 2ll)};
                out_tt.col(0ll) = grid_z.dual_nodes.tail(size_z);
                out_tt.col(1ll) = grid_z.mesh_nodes;
                for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                {
                    out_tt.col(t + 2ll) = states[t].cur_state.col(2ll);
                }
                f << out_tt.format(commaFmt) << '\n';
                f.close();
            }

            // print grids
            {
                ofstream f{std::string{"output/z_grid.csv"}};
                f << grid_z.mesh_nodes.format(commaFmt) << '\n';
                f.close();
            }
            {
                ofstream f{std::string{"output/r_grid.csv"}};
                f << grid_r.mesh_nodes.format(commaFmt) << '\n';
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
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }
#pragma endregion
    }
};

Wrapper::Wrapper(const json &data)
{
    const std::string control_type{data["history"].value<std::string>("control_type", "RFP")};
    if (control_type == "RFP")
    {
        std::cout << "Reservoir flow profile is used as well control condition...\n";
        using Well_t = GPN::Wells::ResFlowProfileControl::WellReservoirFlowProfileControl<History, PhasePropertiesJT, Grids::CylinderGridRock, CrossFlows>;
        WrapperFactory<Well_t>::choose_well(data);
    }
    else if (control_type == "bottomhole_pressure")
    {
        std::cout << "Bottomhole pressure is used as well control condition...\n";
        using Well_t = GPN::Wells::BotHolePresControl::WellBottomHolePressureControl     <History, PhasePropertiesJT, Grids::CylinderGridRock, CrossFlows>;
        WrapperFactory<Well_t>::choose_well(data);
    }
    else if (control_type == "bottomhole_rate")
    {
        std::cout << "Bottomhole rate is used as well control condition...\n";
        using Well_t = GPN::Wells::BotHoleRateControl::WellBottomHoleRateControl         <History, PhasePropertiesJT, Grids::CylinderGridRock, CrossFlows>;
        WrapperFactory<Well_t>::choose_well(data);
    }
    else
        throw std::runtime_error("Incorrect well control type.");
}

std::vector<RealType> Wrapper::get_times() const
{
    return time;
}

Wrapper::~Wrapper()
{
}