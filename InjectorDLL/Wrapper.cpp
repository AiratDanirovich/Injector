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
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/History/History.hpp>
#include <Injector/History/RatesFactory.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Hydrodynamic/Incompressible/IncompressibleFluid.hpp>
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

/// @brief Initial temperature is assumed to be constant
struct FunctorIC : public InitialConditions::ICFunctorBase
{
    FunctorIC(const Logs::Geotherma &geotherma)
        : geotherma{geotherma}
    {
    }
    RealType operator()(const ptrdiff_t z_id, const ptrdiff_t, RealType) const override
    {
        return geotherma(z_id);
    }

protected:
    const Logs::Geotherma &geotherma;
};

template <typename Grid_t_ptr>
auto ICFactory(RealType t_start, const Grid_t_ptr grid, const Logs::Geotherma &geotherma)
{
    return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{geotherma}, t_start)};
}

template <typename Well_t, typename Hydro_t>
struct FunctorBC : public GPN::BoundaryConditions::BCFunctorBase
{
    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
    using ConvectionFieldFactory_t =
        GPN::FaceProperties::IncompressibleRatesFactory<
            Grid2D_t, Well_t, PhasePropertiesJT, Hydro_t>;

    FunctorBC(
        const Logs::IsPermeable &is_permeable,
        const ConvectionFieldFactory_t &flow_field, // volumetric flow rate
        const cptr<const Grid2D_t> grid_ptr)
        : flow_field{flow_field},
          is_permeable{is_permeable},
          grid_ptr{grid_ptr}
    {
    }

    RealType operator()(ptrdiff_t z_id, RealType r, RealType t) const override
    {
        if (r == grid_ptr->second_coord.dual_front())
            return flow_field.get_heat_flow_in_axes2()(z_id, 0ll) * flow_field.get_temperature();

        if (r == grid_ptr->second_coord.dual_back())
            return 0.0;

        assert(false);
        return 0.0;
    }

    RealType operator()(RealType z, ptrdiff_t r_id, RealType t) const override
    {
        if (z == grid_ptr->first_coord.dual_front())
            return flow_field.get_heat_flow_in_axes1()(0ll, r_id) * flow_field.get_temperature();

        if (z == grid_ptr->first_coord.dual_back())
            return 0.0;

        assert(false);
        return 0.0;
    }

protected:
    const ConvectionFieldFactory_t &flow_field;
    const Logs::IsPermeable &is_permeable;
    const cptr<const Grid2D_t> grid_ptr;
};

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
    const RealType t_start,      // start time in seconds
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
    const auto &grid_r{grid2D->second_coord};
    const auto &grid_z{grid2D->first_coord};

    const ExtrudedCasing extr_completion{
        VarExtrudedCasingFactory::create(completion)};
    // collector
    const Logs::Rocks::CoreSampleLogs core_data{
        is_permeable_stencils,
        is_perforated_stencils,
        porosity_stencils,
        permeability_stencils,
        grid_z};
    // make fluid
    const PhasePropertiesJT water{
        FluidFactory::create_water_JT(
            Viscosity{viscosity},
            GPN::Density{density},
            GPN::SpecificHeatCapacity{capacity},
            GPN::HeatConductivity{heat_conductivity_fluid},
            JouleThomson{joule_thomson})};

    const auto RFP_weights{
        Logs::RFPFactory::create_from_container(
            RFP_weights_stencils,
            core_data.is_permeable)};
    const CrossFlows cross_flows{
        RFP_weights, from_coords, to_layers};
    const auto WFP_weights{
        create_WFP(
            core_data.is_perforated,
            RFP_weights,
            cross_flows)};

    const Well_CrossFlow well{RFP_weights, WFP_weights, cross_flows};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        solid_heatconductivity_stencils,
        porosity_stencils,
        water,
        grid2D->first_coord};

    Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};
    heat_props.apply_well(extr_completion, well);

    FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};
    heat_face_props.apply_well(extr_completion, well);

    // history
    const History history{
        HistoryFactory::createFixedRate(time_intervals, well_rates, inlet_temperatures)};
    // external pressure log
    const auto external_pressure{
        Logs::ExtPressureFactory::create(
            ext_pressure_stencils,
            is_permeable_stencils,
            grid_z)};
    // fluid model for the pressure field
    using IncompressibleFluidField_t =
        decltype(IncompressibleFluidField{
            t_start,
            water,
            core_data.permeability,
            external_pressure,
            well,
            grid2D});

    auto ptr_pressure_field{
        make_shared<IncompressibleFluidField_t>(
            t_start,
            water,
            core_data.permeability,
            external_pressure,
            well,
            grid2D)};

    // rates field factory
    FaceProperties::IncompressibleRatesFactory rates_factory{
        ptr_pressure_field,
        grid2D, well, history, water};
    // initial condition
    Logs::Geotherma geotherma{
        Logs::GeothermaFactory::create(
            geotherma_nodes,
            geotherma_vals,
            z_top,
            grid2D->first_coord)};
    const auto initial_state{ICFactory(t_start, grid2D, geotherma)};
    // boundary conditions
    const GPN::BoundaryConditions::BoundaryConditions bc{
        *grid2D,
        std::make_shared<FunctorBC<
            Well_CrossFlow,
            IncompressibleFluidField_t>>(
            core_data.is_permeable, rates_factory, grid2D),
        BoundaryConditions::BoundaryCondition::second};
    // solver
    using Solver_t = decltype(Solver{
        heat_face_props.medium_heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        rates_factory, initial_state,
        bc, t_start});

    auto solver_ptr = std::make_shared<Solver_t>(
        heat_face_props.medium_heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        rates_factory, initial_state,
        bc, t_start);

    const auto &solver{*solver_ptr};

    SolverManager solver_manager{history, solver_ptr};

    solver_manager.run(t_minor_step);

    const auto &[times, states] = solver.solution();
    this->time = times;

    if (!fs::is_directory("output") || !fs::exists("output")) // Check if src folder exists
    {
        fs::create_directory("output"); // create src folder
    }

    try
    {
        ifstream f("separators.json");
        json data = json::parse(f);

        const std::string sep = data["coeff_sep"];

        auto grid{grid2D->second_coord.mesh_nodes};
        grid(0ll) = grid_r.dual_nodes(1ll);

        const Eigen::IOFormat commaFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, sep, sep, "", "", "", "");
        for (auto z{0ll}, layer_id{0ll}; z < grid2D->first_coord.mesh_nodes.size(); ++z)
        {
            if (core_data.is_permeable(z) == 1.0)
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