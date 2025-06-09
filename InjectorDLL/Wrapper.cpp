#include <vector>
#include <cmath>
#include <memory>
#include <numbers>
#include <fstream>
#include <iostream>
#include <filesystem>

#include <InjectorDLL/Wrapper.h>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/History/History.hpp>
#include <Injector/History/RatesFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include <Eigen/Core>

using namespace std;

using namespace GPN;
using namespace GPN::Phases;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

namespace fs = std::filesystem;

/// @brief Initial temperature is assumed to be constant
struct FunctorIC
{
    FunctorIC(const RealType val)
        : val{val}
    {
    }
    RealType operator()(RealType z, RealType r, RealType t_start) const
    {
        return val;
    }

protected:
    const RealType val;
};

template <typename Grid_t_ptr>
auto ICFactory(RealType t_start, const Grid_t_ptr grid, const RealType val)
{
    return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{val}, t_start)};
}

struct FunctorBC : public GPN::BoundaryConditions::BCFunctorBase
{
    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
    using ConvectionFieldFactory_t =
        GPN::FaceProperties::RatesFactory<
            Grid2D_t, Well_KH, PhaseProperties>;

    FunctorBC(
        RealType inlet_temp,
        const Logs::IsPermeable &is_permeable,
        const ConvectionFieldFactory_t &flow_field, // volumetric flow rate
        const cptr<const Grid2D_t> grid_ptr)
        : inlet_temp{inlet_temp},
          flow_field{flow_field},
          is_permeable{is_permeable},
          grid_ptr{grid_ptr}
    {
    }

    RealType operator()(ptrdiff_t z_id, RealType r, RealType t) const override
    {
        if (r == grid_ptr->second_coord.dual_front())
            return flow_field.get_flow_in_axes2()(z_id, 0ll) * flow_field.get_temperature();

        if (r == grid_ptr->second_coord.dual_back())
            return 0.0;

        assert(false);
        return 0.0;
    }

    RealType operator()(RealType z, ptrdiff_t r_id, RealType t) const override
    {
        if (z == grid_ptr->first_coord.dual_front())
            return flow_field.get_flow_in_axes1()(0ll, r_id) * flow_field.get_temperature();

        if (z == grid_ptr->first_coord.dual_back())
            return 0.0;

        assert(false);
        return 0.0;
    }

protected:
    RealType inlet_temp;
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
    // grid
    const RealType rMin,         // m /* typically would be zero */
    const RealType rMax,         // m
    const RealType q,            // --, q >= 1.0 /* step increment factor */
    const RealType r_max_step,   // m /* maximum allowed step in radial direction */
    const RealType zTop,         // m, /* typically would be zero */
    const RealType z_minor_step, // m, /*maximum step within impermeable layers*/
    // eight +1 vectors of the same size
    // values are in SI
    const VR &thickness,                   // meter
    const VR &heatconductivity_stencils,   // Watt/(m*K)
    const VR &porosity,                    // 0.0 < porosity <= 1.0, --
    const VR &permeability_stencils,       // m^2
    const VR &is_permeable,                // {0, 1}, --
    const VR &is_perforated,               // {0, 1}, --
    const VR &solid_density,               // kg/(m^3)
    const VR &solid_specific_heatcapacity, // J/(kg*K)
    const RealType initial_temperature,    // K // should be log in the future
    // temporal grid
    const RealType t_start,          // start time in seconds
    const VR &time_intervals,        // intervals of const rates)
    const RealType t_minor_step,     // time step used for numerical integration
    // well
    const RealType tube_radius,      // m
    const RealType sandface_radius,  // m
    const RealType well_rate,        // ~1.1E-3 m^3/s
    const RealType inlet_temperature // K
)
{
    // adapt stl container to Eigne conteiner
    LogValuesContainer is_permeable_stencils(is_permeable.size());
    std::copy(is_permeable.cbegin(), is_permeable.cend(), is_permeable_stencils.begin());
    LogValuesContainer is_perforated_stencils(is_perforated.size());
    std::copy(is_perforated.cbegin(), is_perforated.cend(), is_perforated_stencils.begin());
    LogValuesContainer solid_density_stencils(solid_density.size());
    std::copy(solid_density.begin(), solid_density.end(), solid_density_stencils.begin());
    LogValuesContainer solid_specific_heatcapacity_stencils(solid_specific_heatcapacity.size());
    std::copy(solid_specific_heatcapacity.begin(), solid_specific_heatcapacity.end(), solid_specific_heatcapacity_stencils.begin());
    LogValuesContainer porosity_stencils(porosity.size());
    std::copy(porosity.begin(), porosity.end(), porosity_stencils.begin());

    // r_stencils
    GPN::WellHoles well_holes{tube_radius, sandface_radius};
    VR r_stencils = well_holes.generate_log_radial_grid(
        rMin, rMax, q, r_max_step);
    // z-refiner
    GPN::Grids::RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
    // the grid itself
    const auto grid2D{
        Grids::CylinderGridFactory::create(
            refiner,
            Grids::Factory::generate_dual_grid_stencils_from_steps(
                zTop, thickness),
            r_stencils)};
    const auto &grid{grid2D->first_coord};
    // collector
    const Logs::Rocks::CoreSampleLogs core_data{
        is_permeable_stencils,
        is_perforated_stencils,
        porosity_stencils,
        permeability_stencils,
        grid};
    // make fluid
    const PhaseProperties water{
        FluidFactory::create_water(
            Viscosity{viscosity},
            Density{density},
            SpecificHeatCapacity{capacity},
            HeatConductivity{heat_conductivity_fluid})};

    const Well_KH well{
        water, core_data.is_permeable, core_data.is_perforated, core_data.permeability};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        heatconductivity_stencils,
        porosity_stencils,
        water,
        grid2D->first_coord};
    Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};
    heat_props.apply_well(well, water);

    const FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};
    // history
    const std::vector<RealType> rates(time_intervals.size(), well_rate);
    const std::vector<RealType> inlet_temperature_set(time_intervals.size(), inlet_temperature);
    const History history{
        HistoryFactory::create(time_intervals, rates, inlet_temperature_set)};
    // rates field factory
    FaceProperties::RatesFactory rates_factory{
        grid2D, well, history, water};
    // initial condition
    const auto initial_state{ICFactory(t_start, grid2D, initial_temperature)};
    // boundary conditions
    const GPN::BoundaryConditions::BoundaryConditions bc{
        *grid2D,
        std::make_shared<FunctorBC>(
            inlet_temperature, core_data.is_permeable, rates_factory, grid2D),
        BoundaryConditions::BoundaryCondition::second};
    // solver
    using Solver_t = decltype(Solver{
        heat_face_props.heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        rates_factory, initial_state,
        bc, t_start});

    auto solver_ptr = std::make_shared<Solver_t>(
        heat_face_props.heat_conductivity,
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

    const std::string sep{", "};
    const Eigen::IOFormat commaFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, sep, sep, "", "", "", "");
    for (auto z{0ll}, layer_id{0ll}; z < grid2D->first_coord.mesh_nodes.size(); ++z)
    {
        if (core_data.is_permeable(z) == 1.0)
        {
            ofstream f{std::string{"output/layer_"} + std::to_string(layer_id) + std::string{".csv"}};

            f << sep << sep << grid2D->second_coord.mesh_nodes.transpose().format(commaFmt) << '\n';
            for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                f << t << sep << times[t] << sep << states[t].cur_state.row(z).format(commaFmt) << '\n';
            }

            f.close();
            ++layer_id;
        }
    }

    {
        ofstream f{std::string{"output/well_temperature.csv"}};
        f << sep << sep << grid2D->first_coord.mesh_nodes.transpose().format(commaFmt) << '\n';
        for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
        {
            f << t << sep << times[t] << sep << states[t].cur_state.col(0ll).format(commaFmt) << '\n';
        }
        f.close();
    }

    {
        ofstream f{std::string{"output/cement_temperature.csv"}};
        f << sep << sep << grid2D->first_coord.mesh_nodes.transpose().format(commaFmt) << '\n';
        for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
        {
            f << t << sep << times[t] << sep << states[t].cur_state.col(1ll).format(commaFmt) << '\n';
        }
        f.close();
    }

    // print grids
    {
        ofstream f{std::string{"output/z_grid.csv"}};
        f << grid2D->first_coord.mesh_nodes.transpose().format(commaFmt) << '\n';
        f.close();
    }
    {
        ofstream f{std::string{"output/r_grid.csv"}};
        f << grid2D->second_coord.mesh_nodes.transpose().format(commaFmt) << '\n';
        f.close();
    }

    {
        ofstream f{std::string{"output/data.txt"}};
        f << "ghost layer height:   " << grid.mesh_nodes(well.ghost_layer_cell_id) << " m" << endl;
        f << "top collector height: " << grid.mesh_nodes(well.top_collector_cell_id) << " m" << endl;
        f.close();
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