#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Properties/FaceProperties.hpp>
#include <Injector/Solver/SplittingMethod/BaseSplit.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::EqSolver::SplittingMethod;

// input data
const auto z_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
const auto r_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};
const auto conductivity_stencils{
    Logs::RawDataFactory::generate_conductivity_StepProperty(z_stencils)};

// Tests Cylinder grid, (r; z)
TEST_CASE("BaseSplitTest")
{
    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
    const auto &grid{grid2D->first_coord};

    const Logs::HeatLogsFactory heat_factory{
        conductivity_stencils, grid};

    Properties::HeatConductivity conductivity_field{
        Properties::FieldFactory::create(
            heat_factory.conductivity,
            grid2D)};

    FaceProperties::HeatConductivity conductivity_field_face{
        FaceProperties::FaceInterpolatedFieldFactory::create(
            conductivity_field,
            grid2D)};

    BaseSplit base_split1{
        conductivity_field_face.face_vals_axes2,
        grid2D->first_coord.mesh_size(),
        grid2D->second_coord.mesh_size()};

    BaseSplit base_split2{
        conductivity_field_face.face_vals_axes1,
        grid2D->second_coord.mesh_size(),
        grid2D->first_coord.mesh_size()};
}
