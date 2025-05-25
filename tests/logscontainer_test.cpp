
#include <Injector/Grids/ConcreteGrids.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Properties/LogsContainer.hpp>

#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::Logs;
using namespace GPN::CoordinateTypes;

TEST_CASE("LogsContainerTest")
{
    const auto grid_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
    const auto is_permeable_stencils{Logs::Factory::generate_is_permeable_StepProperty(grid_stencils)};

    const auto permeability_stencils{Logs::Factory::generate_permeability(grid_stencils, is_permeable_stencils)};
    const auto porosity_stencils{Logs::Factory::generate_porosity_StepProperty(grid_stencils, is_permeable_stencils)};

    const auto grid{ZGrid{GridDual{grid_stencils}}};
    const auto water{FluidFactory::create_water(0.0, 0.0)};

    // Logs::GeologyLogsContainer logs_container{
    //     grid,
    //     water,
    //     is_permeable_stencils,
    //     permeability_stencils,
    //     porosity_stencils,

    // };

    auto is_permeable{
        IsPermeable{
            StepPropertyGrid{
                StepProperty{is_permeable_stencils},
                grid}}};

    const auto permeability{
        Logs::Permeability{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    permeability_stencils},
                grid} * // guarantee that porosity is zero in rocks
                is_permeable,
            is_permeable}};

    const auto porosity{
        Logs::Porosity{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    porosity_stencils},
                grid} * // guarantee that porosity is zero in rocks
                is_permeable,
            is_permeable}};
}
