#include <Injector/Grids/History.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;

// Tests Cylinder grid, (r; z)
TEST_CASE("HistoryTest")
{
    auto time_moments{ 
        Grids::Factory::generate_dual_grid_stencils_uniform(
            0, 90, 91 // linspace operator
        )
    };

    auto rates{
        Logs::Factory::generate_rates_StepProperty(
            time_moments
        )
    };

    auto history{
        History{
            StepProperty{rates}, 
            TemporalGrid{time_moments}
        }
    };

    const auto& grid{history.grid};
    for(auto id{0ll}; id < grid.dual_nodes.size(); ++id)
        CHECK(grid.dual_nodes(id) == grid.dual_stencils(id));
}
