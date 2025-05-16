
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/ConcreteGrids.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Model/Well.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::Phases;
using namespace GPN::CoordinateTypes;

TEST_CASE("Well_KH_Test")
{
    std::vector<RealType> grid_stencils{0.0, 1.0, 3.0, 7.0, 10.0};

    std::vector<RealType> permeability_stencils(grid_stencils.size() - 1ull, 1.0);

    std::vector<RealType> is_permeable_stencils(grid_stencils.size() - 1ull, 1.0);

    auto grid{ZGrid{GridDual{grid_stencils}}};
    auto is_permeable{
        IsPermeable{
            StepPropertyGrid{
                StepProperty{
                    is_permeable_stencils},
                grid}}};

    const auto permeability{
        Logs::Permeability{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    permeability_stencils},
                grid} * // guarantee that porosity is zero in rocks
                is_permeable,
            is_permeable}};

    const Well_KH well{
        FluidFactory::create_water(1.0, 1.0),
        is_permeable, permeability};

    const RealType rate{1.0};
    auto rfp = well.get_RFP(rate);

    for (auto i{0ll}; i < rfp.size(); ++i)
    {
        CHECK(
            rfp(i) ==
            rate / (grid_stencils.back() - grid_stencils.front()) *
                grid.dual_steps(i));
    }
}
