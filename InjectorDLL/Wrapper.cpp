#include <vector>
#include <cmath>
#include <memory>
#include <cassert>

#include <InjectorDLL/Wrapper.h>

#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/FieldsFactory.hpp>
#include <Injector/Properties/FlowField.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

using namespace GPN;
using namespace GPN::Phases;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

struct ExactSolution
{
    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
    ExactSolution(
        const Properties::MediumHeatVolumetricCapacity<Grid2D_t> &volumetric_capacity,
        const Properties::HeatConductivity<Grid2D_t> &heat_conductivity,
        RealType q)
        : kappa{volumetric_capacity, heat_conductivity},
          heat_conductivity{heat_conductivity},
          grid{*heat_conductivity.grid},
          q{q}
    {
    }

    RealType operator()(RealType z, RealType r, RealType t) const
    {
        return -q * heat_conductivity.value(0ull, 0ull) *
               std::expint(-r * r / (4.0 * kappa.value(0ull, 0ull) * t));
    }

    template <typename Grid_t>
    RealType operator()(ptrdiff_t z, ptrdiff_t r, RealType t, const Grid_t &grid) const
    {
        const auto [zv, rv]{grid.coordinates(z, r)};
        return (*this)(zv, rv, t);
    }

protected:
    const Properties::ThermalDiffusivity<Grid2D_t> kappa;
    const Properties::HeatConductivity<Grid2D_t> &heat_conductivity;
    const RealType q;
    const Grid2D_t &grid;
};

struct FunctorIC
{
    FunctorIC(const ExactSolution &es)
        : es{es}
    {
    }
    RealType operator()(RealType z, RealType r, RealType t0) const
    {
        return es(z, r, t0);
    }

protected:
    const ExactSolution &es;
};

template <typename Grid_t_ptr>
auto initialcondition_factory(RealType t0, const Grid_t_ptr grid, const ExactSolution &es)
{
    return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{es}, t0)};
}

struct AFunctorBC : public GPN::BoundaryConditions::BCFunctorBase
{
    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
    AFunctorBC(
        const ExactSolution &es,
        const cptr<Grid2D_t> grid2D)
        : es{es},
          grid2D{grid2D}
    {
    }

    RealType operator()(ptrdiff_t z_id, RealType r, RealType t) const override
    {
        RealType z{grid2D->first_coord.mesh_nodes(z_id)};
        if (r == grid2D->second_coord.dual_front())
            r = grid2D->second_coord.mesh_front();
        else if (r == grid2D->second_coord.dual_back())
            r = grid2D->second_coord.mesh_back();
        else
            assert(false);

        return es(z, r, t);
    }

    RealType operator()(RealType z, ptrdiff_t r_id, RealType t) const override
    {
        RealType r{grid2D->second_coord.mesh_nodes(r_id)};
        if (z == grid2D->first_coord.dual_front())
            z = grid2D->first_coord.mesh_front();
        else if (z == grid2D->first_coord.dual_back())
            z = grid2D->first_coord.mesh_back();
        else
            assert(false);

        return es(z, r, t);
    }

protected:
    const ExactSolution &es;
    const cptr<Grid2D_t> grid2D;
};

Wrapper::Wrapper(
    const RealType q, // = 1.0 heat rate
    // fluid params in SI
    const RealType density,   // kg/(m^3)
    const RealType capacity,  // J/(kg*K) /* specific heat capacity */
    const RealType viscosity, // Pa*s
    // grid
    const RealType rMin,       // m /* typically would be zero */
    const RealType rMax,       // m
    const size_t r_nodes_nmbr, // -- /* number of nodes in r-direction, including first and last ones */
    const RealType zTop,       // m, /* typically would be zero */
    //    const size_t nZ, // = thickness.size()
    // six vectors of the same size
    // values are in SI
    const VR &thickness,                   // meter
    const VR &heatconductivity_stencils,   // Watt/(m*K)
    const VR &porosity,                    // 0.0 < porosity <= 1.0, --
    const VR &is_permeable_stencils,       // {0, 1}, --
    const VR &solid_density,               // kg/(m^3)
    const VR &solid_specific_heatcapacity, // J/(kg*K)
    // temporal grid
    const RealType t_start,  // start time in seconds
                             //    const size_t nt, // = time_intervals.size()
    const VR &time_intervals // intervals of const rates)
)
{
    // adapt stl container to Eigne conteiner
    LogValuesContainer solid_density_stencils(solid_density.size());
    std::copy(solid_density.begin(), solid_density.end(), solid_density_stencils.begin());
    LogValuesContainer solid_specific_heatcapacity_stencils(solid_specific_heatcapacity.size());
    std::copy(solid_specific_heatcapacity.begin(), solid_specific_heatcapacity.end(), solid_specific_heatcapacity_stencils.begin());
    LogValuesContainer porosity_stencils(porosity.size());
    std::copy(porosity.begin(), porosity.end(), porosity_stencils.begin());

    // make grid1D
    const VR z_stencils(
        Grids::Factory::generate_dual_grid_stencils_from_steps(
            zTop, thickness));
    const RealType zBottom{z_stencils.back()};
    const VR r_stencils{
        Grids::Factory::generate_dual_grid_stencils_uniform(
            Segment{rMin, rMax}, r_nodes_nmbr)};
    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
    const auto &grid{grid2D->first_coord};
    // make fluid
    const Water water{
        FluidFactory::create_water(
            Viscosity{viscosity},
            Density{density},
            SpecificHeatCapacity{capacity})};
    // time moments
    const VR t_stencils(
        Grids::Factory::generate_dual_grid_stencils_from_steps(
            t_start, time_intervals));
    // solver
    const Logs::Rocks::IsPermeableLog core_data{
        is_permeable_stencils,
        grid};
    const auto flow_field{
        FaceProperties::FlowFactory::zero_flow(
            core_data.is_permeable, *grid2D)};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        heatconductivity_stencils,
        porosity_stencils,
        Phases::FluidFactory::create_water(1.0, 1.0),
        grid};
    const Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};

    const FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};

    // exact solution
    ExactSolution es{heat_props.medium_vol_heatcapacity,
                     heat_props.heat_conductivity, q};
    // initial conditions
    const auto initial_state{initialcondition_factory(t_start, grid2D, es)};
    // boundary conditions
    const GPN::BoundaryConditions::BoundaryConditions bc{
        *grid2D, std::make_shared<AFunctorBC>(es, grid2D)};

    Solver solver{
        heat_face_props.heat_conductivity,
        flow_field, grid2D,
        heat_props.medium_vol_heatcapacity,
        initial_state,
        bc, t_start};

    // advance in time
    for (size_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
    {
        solver.advance(time_intervals[t_step]);

        const auto &[t, sol] = solver.solution().back();
        // save time
        time.push_back(t);
        // save T
        t_radial_distribution.push_back(std::vector<RealType>(sol.cols(), -1001.0));
        // copy T(r) at z = (zTop + zBottom)/2.0
        for (auto id{0ll}; id < sol.cols(); ++id)
            t_radial_distribution.back()[id] = sol(sol.rows() / 2, id);
    }
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