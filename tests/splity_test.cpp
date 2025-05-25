#include <iostream>
#include <memory>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Properties/FaceProperties.hpp>
#include <Injector/Solver/SplittingMethod/SplitY.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::EqSolver;
using namespace GPN::Grids;
using namespace GPN::EqSolver::SplittingMethod;

using VR = std::vector<RealType>;

// input data
const std::ptrdiff_t nLayers{5ull};
const auto z_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};
const auto r_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};
    
const VR is_permeable_stencils(nLayers, 1.0);
const auto porosity_stencils{
    Logs::RawDataFactory::generate_porosity(z_stencils, is_permeable_stencils)};

// heat logs
const auto solid_density_stencils{
    Logs::RawDataFactory::generate_solid_density(z_stencils)};
const auto solid_specific_heatcapacity_stencils{
    Logs::RawDataFactory::generate_solid_specific_heatcapacity(z_stencils)};
const auto heatconductivity_stencils{
    Logs::RawDataFactory::generate_conductivity(z_stencils)};

TEST_CASE("Solver", "splitY")
{
    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
    const auto &grid{grid2D->first_coord};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        heatconductivity_stencils,
        porosity_stencils,
        Phases::FluidFactory::create_water(1.0, 1.0),
        grid};
    
    const Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};

    const FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};

    const double tol = 1E-14;
    SplitY splity{heat_face_props.heat_conductivity, grid2D};

    for (const auto &m : splity.LaplaceTerms())
    {
        REQUIRE(m.rows() == m.cols());
        {
            auto col{0ll};
            INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << 0.0 << ", r: " << m.coeff(col, col + 1ll));
            CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col + 1ll), tol));
        }
        for (auto col{1ll}; col < m.cols() - 1; ++col)
        {
            INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << m.coeff(col, col - 1ll) << ", r: " << m.coeff(col, col + 1ll));
            CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll) + m.coeff(col, col + 1ll), tol));
            INFO("" << "col: " << col << ", c: " << m.coeff(col, col));
            CHECK(m.coeff(col, col) > tol);
        }
        {
            auto col{m.cols() - 1ll};
            INFO("" << "col: " << col << ", c: " << m.coeff(col, col) << ", l: " << m.coeff(col, col - 1ll) << ", r: " << 0.0);
            CHECK_THAT(-m.coeff(col, col), WithinRel(m.coeff(col, col - 1ll), tol));
        }
    }
}