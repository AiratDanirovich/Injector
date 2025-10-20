#include <memory>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>

#include "includes/BCFunctor.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;

auto z_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
auto r_stencils{Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};

// Tests Cartesian grid
TEST_CASE("PhaseProperties", "Water")
{
  const auto grid2D{CylinderGridFactory::create(z_stencils, r_stencils)};
  BoundaryConditions::GeneralBC bc{grid2D, std::make_shared<BCFunctor>(BCFunctor{1.0})};
}
