#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/FlowField.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::Phases;
using namespace GPN::CoordinateTypes;

std::vector<RealType> z_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
std::vector<RealType> r_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
std::vector<RealType> permeability_stencils(z_stencils.size() - 1ull, 1.0);
std::vector<RealType> is_permeable_stencils(z_stencils.size() - 1ull, 1.0);

TEST_CASE("RFP_reservoir")
{
    const auto grid2D{ Grids::CylinderGridFactory::create(
        z_stencils, r_stencils)};

    const auto is_permeable{
        Logs::IsPermeableFactory::create(
            is_permeable_stencils,
            grid2D->first_coord)};

    const auto permeability{
        Logs::PermeabilityFactory::create(
            permeability_stencils,
            is_permeable_stencils,
            grid2D->first_coord)};

    const auto water{FluidFactory::create_water(1.0, 1.0)};

    const Well_KH well{
        water, is_permeable, permeability};

    const RealType well_rate{1.0};
    const auto rfp{
        RFPFactory::create(
            well.get_RFP(well_rate), is_permeable)};

    for (auto i{0ll}; i < rfp.size(); ++i)
    {
        CHECK(
            rfp(i) ==
            well_rate / (z_stencils.back() - z_stencils.front()) *
                grid2D->first_coord.dual_steps(i));
    }

    CHECK(rfp.log_vals.sum() == well_rate);

    Properties::ReservoirFlowField flow_field{
        rfp, well_rate,
        *grid2D};

    // check the first column of verticle flow
    REQUIRE(flow_field.axes1_as_face_normal.rows() == grid2D->first_coord.dual_size());
    REQUIRE(flow_field.axes1_as_face_normal.cols() == grid2D->second_coord.mesh_size());
    for (auto row{0ll}, col{0ll}; row < flow_field.axes1_as_face_normal.rows(); ++row)
    {
        CHECK(flow_field.axes1_as_face_normal(row, col) == flow_field.axes1_as_face_normal(0, col));
        CHECK(flow_field.axes1_as_face_normal(row, col) == well_rate);
    }

    // check other columns of verticle flow
    for (auto col{1ll}; col < flow_field.axes1_as_face_normal.cols(); ++col)
        for (auto row{0ll}; row < flow_field.axes1_as_face_normal.rows(); ++row)
        {
            CHECK(flow_field.axes1_as_face_normal(row, col) == (RealType)0.0);
        }

    // check columns of horizontal flow
    REQUIRE(flow_field.axes2_as_face_normal.rows() == grid2D->first_coord.mesh_size());
    REQUIRE(flow_field.axes2_as_face_normal.cols() == grid2D->second_coord.dual_size());
    for (auto col{1ll}; col < flow_field.axes2_as_face_normal.cols(); ++col)
        for (auto row{0ll}; row < flow_field.axes2_as_face_normal.rows(); ++row)
        {
            CHECK(flow_field.axes2_as_face_normal(row, col) == flow_field.axes2_as_face_normal(row, 0ll));
        }
}
