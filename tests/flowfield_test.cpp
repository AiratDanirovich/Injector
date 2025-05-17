
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/ConcreteGrids.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Properties/FlowField.hpp>

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

    const RealType well_rate{1.0};
    auto rfp = well.get_RFP(well_rate, z_grid);

    for (auto i{0ll}; i < rfp.size(); ++i)
    {
        CHECK(
            rfp(i) ==
            well_rate / (z_grid_stencils.back() - z_grid_stencils.front()) *
                z_grid.dual_steps(i));
    }

    CHECK(rfp.log_vals.sum() == well_rate);

    Properties::ReservoirFlowField flow_field{
        rfp, well_rate,
        *grid2D};

    // check the first column of verticle flow
    for (auto row{0ll}, col{0ll}; row < flow_field.axes1_as_face_normal.rows(); ++row)
    {
        CHECK(flow_field.axes1_as_face_normal(row, col) == flow_field.axes1_as_face_normal(0, col));
        CHECK(flow_field.axes1_as_face_normal(row, col) == well_rate);
    }

    // check other columns of verticle flow
    for (auto row{0ll}; row < flow_field.axes1_as_face_normal.rows(); ++row)
        for (auto col{1ll}; col < flow_field.axes1_as_face_normal.cols(); ++col)
        {
            CHECK(flow_field.axes1_as_face_normal(row, col) == (RealType)0.0);
        }

    // check columns of horizontal flow
    for (auto row{0ll}; row < flow_field.axes2_as_face_normal.rows(); ++row)
        for (auto col{1ll}; col < flow_field.axes2_as_face_normal.cols(); ++col)
        {
            CHECK(flow_field.axes2_as_face_normal(row, col) == flow_field.axes2_as_face_normal(row, 0ll));
        }
}
