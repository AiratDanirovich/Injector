#include <iostream>
#include <memory>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/FaceProperties.hpp>
#include <Injector/Solver/SolverFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

#include <Injector/Solver/BoundaryConditions.hpp>

#include <catch2/catch_test_macros.hpp>

struct ABCFunctor : public GPN::BoundaryConditions::BCFunctorBase
{
    ABCFunctor(RealType val) : val{val} {}

    RealType operator()(
        RealType x, RealType y, RealType t) const override
    {
        return val;
    }

protected:
    RealType val;
};

using namespace std;
using namespace GPN;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

RealType val{1.0};
const auto z_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
const auto r_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};

// hydrodynamic logs
const auto is_permeable_stencils{
    Logs::RawDataFactory::generate_is_permeable(z_stencils)};
const auto porosity_stencils{
    Logs::RawDataFactory::generate_porosity(z_stencils, is_permeable_stencils)};
const auto permeability_stencils{
    Logs::RawDataFactory::generate_permeability(z_stencils, is_permeable_stencils)};
// heat logs
const auto solid_density_stencils{
    Logs::RawDataFactory::generate_solid_density(z_stencils)};
const auto solid_specific_heatcapacity_stencils{
    Logs::RawDataFactory::generate_solid_specific_heatcapacity(z_stencils)};
const auto heatconductivity_stencils{
    Logs::RawDataFactory::generate_conductivity_StepProperty(z_stencils)};

TEST_CASE("Solver")
{
    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
    const auto &grid{grid2D->first_coord};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        heatconductivity_stencils,
        porosity_stencils,
        Phases::FluidFactory::create_water(1.0, 1.0),
        grid};

    const Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};

    FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};

    const Problem::InitialCondition initial_state{
        State::State2D::FillWithConst(
            *grid2D, val)};

    const BoundaryConditions::BoundaryConditions bc{
        *grid2D,
        make_shared<ABCFunctor>(val)};

    const RealType well_rate{1.0};

    const auto flow_field{solver_factory.flow_field(well_rate)};

    Solver solver{
        heat_face_props.heat_conductivity,
        flow_field, grid2D,
        heat_props.medium_vol_heatcapacity,
        initial_state,
        bc, 0.0};

    solver.advance(0.005);
}