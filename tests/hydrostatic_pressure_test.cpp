#include <iostream>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>

#include <Injector/Model/Phases/FluidFactory.hpp>

#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/transfer_to_vector.hpp"
#include "includes/generate_stencils_and_steps.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace std;
using namespace GPN;
using namespace GPN::Logs;
using namespace Catch;
using namespace Catch::Matchers;

/*input data*/
// z-grid data
const auto nLayers{5ull};
const auto z_grid_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, nLayers)};

const RealType z_top{0.0};

const RealType z_ref{1000};

const auto fluid{Phases::FluidFactory::create_water(0.0, 0.0)};

const RealType P_bot{Gravity::value()*fluid.density*z_ref+1};

const RealType slope{Gravity::value()*fluid.density};

TEST_CASE("HydrostaticPressureTest")
{
    {
        const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(z_grid_stencils)};
        const Logs::HydrostaticPressureFactory hydrostatic_factory{z_ref, fluid, grid};
        const auto hydrostatic_pressure{hydrostatic_factory.create(P_bot)};

        const auto &mesh{grid.mesh_nodes};
        CHECK(hydrostatic_pressure.size() == grid.mesh_nodes.size());
        for (auto i{0ll}; i < hydrostatic_pressure.size(); ++i)
            CHECK_THAT(P_bot + slope * (mesh(i) - z_ref), WithinRel(hydrostatic_pressure(i), 1E-10));
    }
    {
        const auto grid{Grids::Factory::create_axes<CoordinateTypes::Z>(z_grid_stencils)};
        const Logs::HydrostaticPressureFactory hydrostatic_factory{z_ref, fluid, grid};
        const auto hydrostatic_pressure{hydrostatic_factory.create_const(P_bot)};

        const auto &mesh{grid.mesh_nodes};
        CHECK(hydrostatic_pressure.size() == grid.mesh_nodes.size());
        for (auto i{0ll}; i < hydrostatic_pressure.size(); ++i)
            CHECK_THAT(P_bot, WithinRel(hydrostatic_pressure(i), 1E-10));
    }
}
