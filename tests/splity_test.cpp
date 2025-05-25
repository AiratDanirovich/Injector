#include <iostream>
#include <memory>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Solver/SplittingMethod/SplitY.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;

using namespace GPN;
using namespace GPN::EqSolver;
using namespace GPN::Grids;
using namespace GPN::EqSolver::SplittingMethod;

// input data
const auto z_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
const auto r_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};
const auto conductivity_stencils{
    Logs::RawDataFactory::generate_conductivity_StepProperty(z_stencils)};

TEST_CASE("Solver", "splitY")
{
    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
    const auto &grid{grid2D->first_coord};
    const Logs::HeatLogsFactory heat_factory{conductivity_stencils, grid};
    const Properties::HeatConductivity conductivity_field{
        Properties::FieldFactory::create(
            heat_factory.conductivity,
            grid2D)};

    const double tol = 1E-14;
    SplitY splity{conductivity_field, grid2D};

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