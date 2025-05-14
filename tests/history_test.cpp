
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
    auto time_moments{
        Grids::Factory::generate_dual_grid_stencils_uniform(
            0, 31536000, 12 // seconds in 12 months
            )};

            std::cout << "fff\n";

    auto time_steps{
        DualStepsContainer{Grids::Factory::generate_dual_grid_steps(
            time_moments)}};

    auto time_grid{
        TemporalGridDual{GridDualStencils{time_steps}}};
#pragma endregion

    auto rates{
        InjectorRate{
            StepPropertyGrid{
                Logs::Factory::generate_rates_StepProperty(
                    time_moments),
                time_grid}}};

    auto history{
        History{rates}};

    // const auto &grid{history.grid};
    // for (auto id{0ll}; id < grid.dual_nodes.size(); ++id)
    //     CHECK(grid.dual_nodes(id) == grid.dual_stencils(id));
}
