#include <iostream>
#include <memory>

#include <catch2/catch_test_macros.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Grids/PhysicalField.hpp>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <Injector/Solver/SplittingMethod/Solver.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::State;
using namespace GPN::EqSolver::Problem;
using namespace GPN::EqSolver::SplittingMethod;

struct BCFunctor : public GPN::BoundaryConditions::BCFunctorBase
{
    RealType operator()(
        RealType x, RealType y, RealType t) const override
    {
        return 1.0;
    }
};

TEST_CASE("Solver")
{
    auto z_grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
    auto z_grid{ZGrid{GridDual{z_grid_stencils}}};
#pragma region GRID_2D
    // generate 1D grids in every direction --- points of property jumps
    StructuredCylinderGrid2DAxisymmetric
        grid2D{z_grid,
               RGrid{
                   GridDual{
                       GridDualStencils{
                           Grids::Factory::
                               generate_dual_grid_stencils_uniform(0, 1, 11)}}}};
#pragma endregion
#pragma region HEAT-CONDUCTIVITY
    // generate heat conductivity field
    Properties::HeatConductivity conductivity_field{
        Logs::HeatConductivity{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    Logs::Factory::generate_conductivity_StepProperty(
                        grid2D.first_coord.dual_stencils)},
                grid2D.first_coord}},
        grid2D};
#pragma endregion

    auto solid_density_stencils{
        Logs::Factory::generate_solid_density_StepProperty(
            z_grid_stencils)};
    auto solid_density{
        SolidDensity{
            StepPropertyGrid{
                StepProperty{
                    solid_density_stencils},
                z_grid}}};

    auto solid_specific_heatcapacity_stencils{
        Logs::Factory::generate_solid_specific_heatcapacity_StepProperty(
            z_grid_stencils)};
    auto solid_specific_heatcapacity{
        SolidSpecificHeatCapacity{
            StepPropertyGrid{
                StepProperty{
                    solid_specific_heatcapacity_stencils},
                z_grid}}};

    auto solid_volumetric_heatcapacity{
        SolidVolumetricHeatCapacity{
            solid_density, solid_specific_heatcapacity}};

    auto water{Phases::FluidFactory::create_water(300, 10)};
    auto is_permeable{
        IsPermeable{
            StepPropertyGrid{
                StepProperty{
                    Logs::Factory::generate_is_permeable_StepProperty(z_grid_stencils)},
                z_grid}}};

    auto porosity_stencils{Logs::Factory::generate_porosity_StepProperty(z_grid_stencils)};
    auto porosity{
        Porosity{
            StepPropertyGrid{
                StepProperty{
                    porosity_stencils} * // guarantee that porosity is zero in rocks
                    is_permeable,
                z_grid},
            is_permeable}};

    auto capacity{
        HeatVolumetricCapacity{
            porosity, solid_volumetric_heatcapacity, water}};

    Properties::HeatVolumetricCapacity capacity_field{
        capacity,
        grid2D};

    InitialCondition
        initial_state{
            State2D::FillWithConst(
                grid2D, 1.0)};

    GPN::BoundaryConditions::BoundaryConditions bc{grid2D, std::make_shared<BCFunctor>(BCFunctor{})};

    const double tol = 1E-8;

    Solver solver{
        conductivity_field, grid2D,
        capacity_field,
        initial_state,
        bc, 0.0};
}