
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/FieldsFactory.hpp>

#include <Injector/Properties/FaceProperties.hpp>

#include <catch2/catch_test_macros.hpp>
#include <iostream>

using namespace std;
using namespace GPN;
using namespace GPN::CoordinateTypes;

/*input data*/
// z-grid data
const auto nLayers{5ull};
const auto grid_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};
// hydrodynamic logs
const auto is_permeable_stencils{
    Logs::RawDataFactory::generate_is_permeable_const(1.0, grid_stencils)};
const auto porosity_stencils{
    Logs::RawDataFactory::generate_porosity(grid_stencils, is_permeable_stencils)};
const auto permeability_stencils{
    Logs::RawDataFactory::generate_permeability(grid_stencils, is_permeable_stencils)};

const auto ext_pressure_stencils{
    Logs::RawDataFactory::generate_ext_pressure(grid_stencils, is_permeable_stencils)};
const auto skin_stencils{
    Logs::RawDataFactory::generate_skin(grid_stencils, is_permeable_stencils)};
// heat logs
const auto solid_density_stencils{
    Logs::RawDataFactory::generate_solid_density(grid_stencils)};
const auto solid_specific_heatcapacity_stencils{
    Logs::RawDataFactory::generate_solid_specific_heatcapacity(grid_stencils)};
const auto heatconductivity_stencils{
    Logs::RawDataFactory::generate_conductivity(grid_stencils)};

TEST_CASE("InterpolatorTest")
{
    // create 2D grid
    const auto grid2D{Grids::CylinderGridFactory::create(grid_stencils, grid_stencils)};
    const auto &grid_z{grid2D->first_coord};
    const auto &grid_r{grid2D->second_coord};

    const Logs::Permeability permeability{
        Logs::PermeabilityFactory::create(
            permeability_stencils,
            is_permeable_stencils,
            grid_z)};

    cout << "permeability:\n";
    std::cout << permeability.log_vals << std::endl
              << std::endl;

    const auto perm_z{Logs::FaceInterpolator::interpolate1D_z(permeability, grid_z)};
    cout << "perm_z:\n";
    cout << perm_z << endl
         << endl;

    for (auto row{0ll}; row < perm_z.rows(); ++row)
    {
        // if (
        //     (std::abs(
        //          grid_z.dual_nodes(row + 1ll) -
        //          (grid_z.mesh_nodes(row) +
        //           grid_z.mesh_nodes(row + 1ll)) /
        //              2.0) < 1E-10)

        //     &&
        //     (std::abs(
        //          grid_z.dual_nodes(0 + 1ll) -
        //          (grid_z.mesh_nodes(0) +
        //           grid_z.mesh_nodes(0 + 1ll)) /
        //              2.0) < 1E-10))
        // {
        //     INFO("" << "row: " << row << ", mesh_nodes: {" << grid_z.mesh_nodes(row) << ", " << grid_z.mesh_nodes(row + 1ll) << "}, dual_node: " << grid_z.dual_nodes(row + 1ll));
        //     CHECK(perm_z(0ll) == perm_z(row));
        // }

        CHECK(perm_z(row) == CartesianCoordinate::face_interpolator(
                                 grid_z.mesh_nodes(row),
                                 grid_z.mesh_nodes(row + 1ll),
                                 grid_z.dual_nodes(row + 1ll),
                                 permeability.log_vals(row),
                                 permeability.log_vals(row + 1ll)));
    }

    const auto perm_r{
        Logs::FaceInterpolator::interpolate1D_r(
            permeability.log_vals.transpose().eval(), grid_r)};
    cout << "perm_r:\n";
    cout << perm_r << endl
         << endl;
    for (auto row{0ll}; row < perm_r.rows(); ++row)
    {
        CHECK(perm_r(row) == R_CylCoord::face_interpolator(
                                 grid_r.mesh_nodes(row),
                                 grid_r.mesh_nodes(row + 1ll),
                                 grid_r.dual_nodes(row + 1ll),
                                 permeability.log_vals(row),
                                 permeability.log_vals(row + 1ll)));
    }

    const Properties::Permeability perm_field{
        Properties::FieldFactory::create(permeability, grid2D)};
    cout << "perm_field:\n";
    cout << perm_field.values() << endl
         << endl;

    REQUIRE(perm_field.values().rows() == permeability.log_vals.rows());
    for (auto row{0ll}; row < perm_field.values().rows(); ++row)
    {
        for (auto col{0ll}; col < perm_field.values().cols(); ++col)
        {
            CHECK(permeability.log_vals(row) == perm_field.value(row, col));
        }
    }

    const auto perm_field_z{
        FaceProperties::FaceInterpolator::interpolate2D_z(
            perm_field.values(), *grid2D)};
    cout << "perm_field_z:\n";
    cout << perm_field_z << endl
         << endl;

    REQUIRE(perm_field_z.rows() == perm_z.rows());
    for (auto row{0ll}; row < perm_field_z.rows(); ++row)
    {
        CHECK(perm_field_z(row, 0ll) == perm_z(row));
        for (auto col{1ll}; col < perm_field_z.cols(); ++col)
        {
            CHECK(perm_field_z(row, 0ll) == perm_field_z(row, col));
        }
    }

    const auto perm_field_r{
        FaceProperties::FaceInterpolator::interpolate2D_r(
            perm_field.values(), *grid2D)};
    cout << "perm_field_r:\n";
    cout << perm_field_r << endl
         << endl;
}
