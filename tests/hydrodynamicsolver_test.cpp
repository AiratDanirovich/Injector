
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::CoordinateTypes;

TEST_CASE("HydrodynamicsSolverTest")
{
#pragma region GRID_2D
    // generate 1D grids in every direction --- points of property jumps
    auto z_grid{
        ZGrid{
            GridDual{
                GridDualStencils{
                    Grids::Factory::
                        generate_dual_grid_stencils_uniform(0, 1, 5)}}}};
    auto r_grid{
        RGrid{
            GridDual{
                GridDualStencils{
                    Grids::Factory::
                        generate_dual_grid_stencils_uniform(0, 1, 11)}}}};
    cptr<StructuredCylinderGrid2DAxisymmetric>
        grid2D{std::make_shared<StructuredCylinderGrid2DAxisymmetric>(
            z_grid, r_grid)};
#pragma endregion
#pragma region IS_PERMEABLE
    const auto is_permeable_stencils{
        Logs::Factory::generate_is_permeable_StepProperty(grid2D->first_coord.dual_stencils)};
    auto is_permeable{
        Logs::IsPermeable{
            Logs::StepPropertyGrid{
                Logs::StepProperty{is_permeable_stencils},
                grid2D->first_coord}}};
#pragma endregion

#pragma region PERMEABILITY
    // generate permeability log
    const auto permeability{
        Logs::Permeability{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    Logs::Factory::generate_permeability_StepProperty(
                        grid2D->first_coord.dual_stencils, is_permeable_stencils)},
                grid2D->first_coord},
            is_permeable}};

    // generate permeability 2D field
    Properties::Permeability permeability_field{
        permeability,
        grid2D};
#pragma endregion
#pragma region POROSITY
    // generate porosity log

    const auto porosity{
        Logs::Porosity{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    Logs::Factory::generate_porosity_StepProperty(
                        grid2D->first_coord.dual_stencils, is_permeable_stencils)},
                grid2D->first_coord},
            is_permeable}};

    // generate porosity 2D field
    Properties::Porosity porosity_field{
        porosity,
        grid2D};
#pragma endregion
#pragma region HEAT-CONDUCTIVITY
    // generate heat conductivity log
    Logs::HeatConductivity conductivity{
        Logs::StepPropertyGrid{
            Logs::StepProperty{
                Logs::Factory::generate_conductivity_StepProperty(
                    grid2D->first_coord.dual_stencils)},
            grid2D->first_coord}};
    // another way to generate logs
    // auto conductivity2{
    //     Logs::generate_log<Logs::HeatConductivity>(
    //         Logs::Factory::generate_conductivity_StepProperty(
    //             grid2D.first_coord.dual_stencils),
    //         grid2D.first_coord)};

    Properties::HeatConductivity conductivity_field{
        conductivity,
        grid2D};
#pragma endregion
}
