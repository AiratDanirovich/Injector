
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::CoordinateTypes;

TEST_CASE("LogsTest")
{
    auto grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
    const auto is_permeable_stencils{Logs::Factory::generate_is_permeable_StepProperty(grid_stencils)};

    auto permeability_stencils{Logs::Factory::generate_permeability_StepProperty(grid_stencils, is_permeable_stencils)};

    auto porosity_stencils{Logs::Factory::generate_porosity_StepProperty(grid_stencils, is_permeable_stencils)};

    auto grid{ZGrid{GridDual{grid_stencils}}};
    auto is_permeable{
        IsPermeable{
            StepPropertyGrid{
                StepProperty{
                    Logs::Factory::generate_is_permeable_StepProperty(grid_stencils)},
                grid}}};

    const auto permeability{
        Logs::Permeability{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    permeability_stencils},
                grid},
            is_permeable}};

    const auto porosity{
        Logs::Porosity{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    porosity_stencils},
                grid},
            is_permeable}};
}
