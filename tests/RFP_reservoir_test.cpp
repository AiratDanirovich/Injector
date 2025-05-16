
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/ConcreteGrids.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Model/Well.hpp>

#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::Phases;
using namespace GPN::CoordinateTypes;

TEST_CASE("RFP_reservoir")
{
#pragma region GRID_2D
    std::vector<RealType> z_grid_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
    std::vector<RealType> r_grid_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
    // generate 1D grids in every direction --- points of property jumps
    auto z_grid{ZGrid{GridDual{z_grid_stencils}}};
    auto r_grid{RGrid{GridDual{r_grid_stencils}}};
    cptr<StructuredCylinderGrid2DAxisymmetric>
        grid2D{std::make_shared<StructuredCylinderGrid2DAxisymmetric>(
            z_grid, r_grid)};
#pragma endregion

    std::vector<RealType> permeability_stencils(z_grid_stencils.size() - 1ull, 1.0);

    std::vector<RealType> is_permeable_stencils(z_grid_stencils.size() - 1ull, 1.0);

    auto is_permeable{
        IsPermeable{
            StepPropertyGrid{
                StepProperty{
                    is_permeable_stencils},
                z_grid}}};

    const auto permeability{
        Logs::Permeability{
            Logs::StepPropertyGrid{
                Logs::StepProperty{
                    permeability_stencils},
                z_grid} * // guarantee that porosity is zero in rocks
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
            rate / (z_grid_stencils.back() - z_grid_stencils.front()) *
                z_grid.dual_steps(i));
    }

    CHECK(rfp.log_vals.sum() == rate);

    Properties::RFP rfp_field{
        rfp,
        grid2D};
}
