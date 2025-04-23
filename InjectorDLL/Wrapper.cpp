
#include <vector>

#include <InjectorDLL/Wrapper.h>

#include <Injector/Grids/Factory.hpp>
#include <Injector/Grids/PhysicalField.hpp>
#include <Injector/Properties/Factory.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>

#include <Injector/Solver/SplittingMethod/Solver.hpp>


using namespace GPN;
using namespace GPN::Phases;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

struct ExactSolution
{
    ExactSolution(
        const Properties::HeatVolumetricCapacity &volumetric_capacity,
        const Properties::HeatConductivity &heat_conductivity,
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
    const Properties::ThermalDiffusivity kappa;
    const Properties::HeatConductivity &heat_conductivity;
    const RealType q;
    const Properties::HeatConductivity::Grid_type &grid;
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

struct FunctorBC : public GPN::BoundaryConditions::BCFunctorBase
{
    FunctorBC(const ExactSolution &es)
        : es{es}
    {
    }
    RealType operator()(RealType z, RealType r, RealType t) const override
    {
        return es(z, r, t);
    }

protected:
    const ExactSolution &es;
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
    const VR &conductivity,                // Watt/(m*K)
    const VR &porosity,                    // --
    const VR &is_permeable,                // {0, 1}, --
    const VR &solid_density,               // kg/(m^3)
    const VR &solid_specific_heatcapacity, // J/(kg*K)
    // temporal grid
    const RealType t_start,  // start time in seconds
                             //    const size_t nt, // = time_intervals.size()
    const VR &time_intervals // intervals of const rates)
)
{
    // make grid1D
    const VR z_stencils(
        Grids::Factory::generate_dual_grid_stencils_from_steps(
            zTop, thickness));
    const RealType zBottom{z_stencils.back()};
    const VR r_stencils{
        Grids::Factory::generate_dual_grid_stencils_uniform(
            Segment{rMin, rMax}, r_nodes_nmbr)};
    // make grid2D
    const auto grid2D{
        Grids::Factory::create_cylinder_grid_2D_ptr(
            z_stencils, r_stencils)};
    // heat conductivity
    const Properties::HeatConductivity conductivity_field{
        Properties::Factory::generate_heatconductivity_Property(
            conductivity, grid2D)};
    // make fluid
    const Water water{
        FluidFactory::create_water(
            Viscosity{viscosity},
            Density{density},
            SpecificHeatCapacity{capacity})};
    // volumetric heat capacity of multiphaase system
    const Properties::HeatVolumetricCapacity capacity_field{
        Properties::Factory::generate_volumetric_heatcapacity_Property(
            is_permeable, porosity,
            solid_density, solid_specific_heatcapacity,
            water, grid2D)};

    // exact solution
    ExactSolution es{capacity_field,
                     conductivity_field, q};
    // initial conditions
    const auto initial_state{initialcondition_factory(t_start, grid2D, es)};
    // boundary conditions
    const GPN::BoundaryConditions::BoundaryConditions bc{
        *grid2D, std::make_shared<FunctorBC>(es)};

    // solver
    Solver solver{
        conductivity_field, grid2D,
        capacity_field,
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