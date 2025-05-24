#include <Injector/Grids/Defines.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::Phases;
using namespace GPN::CoordinateTypes;

std::vector<RealType> grid_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
std::vector<RealType> permeability_stencils(grid_stencils.size() - 1ull, 1.0);
std::vector<RealType> is_permeable_stencils(grid_stencils.size() - 1ull, 1.0);

TEST_CASE("Well_KH_Test")
{
    const auto z_grid{
        Grids::Factory::create_axes<CoordinateTypes::Z>(
            grid_stencils)};

    const auto is_permeable{
        Logs::IsPermeableFactory::create(
            is_permeable_stencils,
            z_grid)};

    const auto permeability{
        Logs::PermeabilityFactory::create(
            permeability_stencils,
            is_permeable_stencils,
            z_grid)};

    const auto water{FluidFactory::create_water(1.0, 1.0)};

    const Well_KH well{
        water, is_permeable, permeability};

    const RealType rate{1.0};
    auto rfp = RFPFactory::create(well.get_RFP(rate), is_permeable);

    for (auto i{0ll}; i < rfp.size(); ++i)
    {
        CHECK(
            rfp(i) ==
            rate / (grid_stencils.back() - grid_stencils.front()) *
                z_grid.dual_steps(i));
    }

    CHECK(rfp.log_vals.sum() == rate);
}
