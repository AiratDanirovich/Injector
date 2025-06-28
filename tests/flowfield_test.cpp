#include <iostream>
#include <limits>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/FlowField.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace std;

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::Phases;
using namespace GPN::CoordinateTypes;

using namespace Catch;
using namespace Catch::Matchers;

using VR = std::vector<RealType>;
vector<RealType> make_steps(const VR &data)
{
    VR out(data.size() - 1ull);
    for (auto i{1ull}; i < data.size(); ++i)
        out[i - 1ull] = data[i] - data[i - 1ull];

    return out;
}

LogValuesContainer transfer_to_eigen(const VR &data, const RealType factor = 1.0)
{
    LogValuesContainer out(data.size());
    for (auto i{0ull}; i < data.size(); ++i)
        out(i) = factor * data[i];
    return out;
}

std::vector<RealType> z_stencils{-2.0, -1.0, 0.0, 1.0, 3.0, 7.0, 10.0};
LogValuesContainer z_thickness{transfer_to_eigen(make_steps(z_stencils))};
std::vector<RealType> r_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
std::vector<RealType> permeability_stencils(z_stencils.size() - 1ull, 1.0);
std::vector<RealType> is_permeable_stencils(z_stencils.size() - 1ull, 0.0);

const RealType well_rate{1.0};

const RealType rMax{300.0}; // m
/*well*/
const RealType sandface_radius{0.3}; // m
const RealType column_radius{0.2}; // m
const RealType tube_radius{0.1};     // m

struct Record
{
    const RealType rate, pressure;
} history_record{well_rate, std::numeric_limits<double>::quiet_NaN()};

TEST_CASE("RFP_reservoir")
{
    is_permeable_stencils[0] = 0.0;
    is_permeable_stencils[2] = 0.0;
    is_permeable_stencils[3] = 0.0;

    
    is_permeable_stencils[5] = 1.0;

    std::vector<RealType> is_perforated_stencils{is_permeable_stencils};
    std::vector<RealType> permeability_stencils{is_permeable_stencils};

    auto it = std::ranges::find_if(
        is_perforated_stencils,
        [](RealType v)
        { return v == 1.0; });
//    (*it) = 0.0;

    const auto grid2D{Grids::CylinderGridFactory::create(
        z_stencils, r_stencils)};

    const auto is_permeable{
        Logs::IsPermeableFactory::create(
            is_permeable_stencils,
            grid2D->first_coord)};

    const auto is_perforated{
        Logs::IsPerforatedFactory::create(
            is_perforated_stencils,
            is_permeable_stencils,
            grid2D->first_coord)};

    const auto permeability{
        Logs::PermeabilityFactory::create(
            permeability_stencils,
            is_permeable_stencils,
            grid2D->first_coord)};

    const auto water{FluidFactory::create_water(1.0, 1.0)};

    WellHoles well_holes{TubeInnerRadius{tube_radius}, SandfaceRadius{sandface_radius}};

    const Well_KH well{
        water, is_permeable, is_perforated, permeability, well_holes, rMax};

    {
        const auto rfp{
            RFPFactory::create_from_container(
                well.get_RFP(history_record), is_permeable)};
                
        cout << "rfp:\n"
             << rfp.log_vals.transpose() << endl;
        cout << "top_collector_cell_id:\n"
             << well.top_collector_cell_id() << endl;
             

        for (auto i{0ll}; i < rfp.size(); ++i)
        {
            CHECK(
                rfp(i) ==
                is_permeable(i) * well_rate / (z_thickness * is_permeable.log_vals).sum() *
                    grid2D->first_coord.dual_steps(i));
        }

        CHECK(rfp.log_vals.sum() == well_rate);

        FaceProperties::ReservoirFlowField flow_field{
            FaceProperties::FlowFactory::create(
                rfp, history_record,
                *grid2D)};

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

    {
        FaceProperties::ReservoirFlowField flow_field{
            FaceProperties::FlowFactory::create_from_well(
                history_record, well,
                *grid2D)};


        cout << "horizontl flow:\n"
             << flow_field.axes2_as_face_normal << endl;
        cout << "vertical flow:\n"
             << flow_field.axes1_as_face_normal << endl;

        for (auto col{0ll}; col < grid2D->second_coord.mesh_size(); ++col)
        {
            for (auto row{0ll}; row < grid2D->first_coord.mesh_size(); ++row)
            {
                CHECK_THAT(
                    flow_field.axes2_as_face_normal(row, col) - flow_field.axes2_as_face_normal(row, col + 1ll),
                    WithinRel(
                        -flow_field.axes1_as_face_normal(row, col) + flow_field.axes1_as_face_normal(row + 1ll, col),
                        1E-11));
            }
        }
        // {
        //     // verify the correct flow in cement
        //     const ptrdiff_t col{1ll};
        //     for (auto row{0ll}; row <= well.ghost_layer_cell_id; ++row)
        //     {
        //         CHECK(flow_field.axes1_as_face_normal(row, col) == 0.0);
        //     }
        //     for (auto row{well.ghost_layer_cell_id + 1ll}; row <= well.top_collector_cell_id; ++row)
        //     {
        //         CHECK(flow_field.axes1_as_face_normal(row, col) == flow_field.axes1_as_face_normal(well.ghost_layer_cell_id + 1ll, col));
        //         CHECK(flow_field.axes1_as_face_normal(row, col) < 0.0);
        //     }
        //     for (auto row{well.top_collector_cell_id + 1ll}; row < flow_field.axes1_as_face_normal.rows(); ++row)
        //     {
        //         CHECK(flow_field.axes1_as_face_normal(row, col) == 0.0);
        //     }
        // }

        {
            // verify the correct verticle flow in rocks
            for (auto col{2ll}; col < flow_field.axes1_as_face_normal.cols(); ++col)
            {
                for (auto row{0ll}; row < flow_field.axes1_as_face_normal.rows(); ++row)
                {
                    CHECK(flow_field.axes1_as_face_normal(row, col) == 0.0);
                }
            }
        }
    }
}
