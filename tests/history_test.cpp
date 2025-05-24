
#include <Injector/History/TemporalGrid.hpp>
#include <Injector/History/History.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <catch2/catch_test_macros.hpp>
#include <iostream>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;

// Tests Cylinder grid, (r; z)
TEST_CASE("HistoryTest")
{
#pragma region MAKE-TIME-GRID
    const auto time_moments{
        Grids::Factory::generate_dual_grid_stencils_uniform(
            0, 31536000, 13 // seconds in 12 months
            )};

    const auto time_steps{
        DualStepsContainer{Grids::Factory::generate_dual_grid_steps(
            time_moments)}};

    const auto time_grid{
        TemporalGridDual{GridDualStencils{time_steps}}};
#pragma endregion

    const auto rates{
        InjectorRate{
            StepPropertyGrid{
                Logs::RawDataFactory::generate_rates(
                    time_moments),
                time_grid}}};

    const auto history{
        History{rates}};

    const auto &grid{history.rates.grid};
    for (auto id{0ll}; id < grid.dual_nodes.size(); ++id)
        CHECK(grid.dual_nodes(id) == grid.dual_stencils(id));
}
