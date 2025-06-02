#include <vector>
#include <cmath>
#include <memory>
#include <numbers>

#include <InjectorDLL/Wrapper.h>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
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

#include <Injector/Properties/FieldsFactory.hpp>

using namespace GPN;
using namespace GPN::Phases;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

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

struct FunctorBC : public BoundaryConditions::BCFunctorBase
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
            return flow_field.get_flow_in_axes2()(z_id, 0ll) * inlet_temp;

        if (r == grid_ptr->second_coord.dual_back())
            return 0.0;

        assert(false);
        return 0.0;
    }

    RealType operator()(RealType z, ptrdiff_t r_id, RealType t) const override
    {
        if (z == grid_ptr->first_coord.dual_front())
            return flow_field.get_flow_in_axes1()(0ll, r_id) * inlet_temp;

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
    const RealType rMin, // m /* typically would be zero */
    const RealType rMax, // m
    const size_t rNodes, // -- /* number of nodes in r-direction, including first and last ones */
    const RealType zTop, // m, /* typically would be zero */
    // seven +1 vectors of the same size
    // values are in SI
    const VR &thickness,                   // meter
    const VR &heatconductivity_stencils,   // Watt/(m*K)
    const VR &porosity,                    // 0.0 < porosity <= 1.0, --
    const VR &permeability_stencils,       // m^2
    const VR &is_permeable_stencils,       // {0, 1}, --
    const VR &solid_density,               // kg/(m^3)
    const VR &solid_specific_heatcapacity, // J/(kg*K)
    const RealType initial_temperature,    // K // should be log in the future
    // temporal grid
    const RealType t_start,      // start time in seconds
    const VR &time_intervals,    // intervals of const rates)
    const RealType t_minor_step, // time step used for numerical integration
    // well
    const RealType well_rate,        // ~1.1E-3 m^3/s
    const RealType inlet_temperature // K
)
{
    // adapt stl container to Eigne conteiner
    LogValuesContainer solid_density_stencils(solid_density.size());
    std::copy(solid_density.begin(), solid_density.end(), solid_density_stencils.begin());
    LogValuesContainer solid_specific_heatcapacity_stencils(solid_specific_heatcapacity.size());
    std::copy(solid_specific_heatcapacity.begin(), solid_specific_heatcapacity.end(), solid_specific_heatcapacity_stencils.begin());
    LogValuesContainer porosity_stencils(porosity.size());
    std::copy(porosity.begin(), porosity.end(), porosity_stencils.begin());

    // make grid2D
    const auto grid2D{
        Grids::CylinderGridFactory::create(
            Grids::Factory::generate_dual_grid_stencils_from_steps(
                zTop, thickness),
            Grids::Factory::generate_dual_grid_stencils_uniform(
                Segment{rMin, rMax}, rNodes))};
    const auto &grid{grid2D->first_coord};
    // collector
    const Logs::Rocks::CoreSampleLogs core_data{
        is_permeable_stencils,
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
        water, core_data.is_permeable, core_data.permeability};

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
    const History history{
        HistoryFactory::create(time_intervals, rates)};
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
}

std::vector<RealType> Wrapper::get_times() const
{
    return time;
}
std::vector<std::vector<RealType>> Wrapper::get_temps() const
{
    return t_radial_distribution;
}

Wrapper::~Wrapper()
{
}