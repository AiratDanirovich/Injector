#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Phases;

// Tests Cartesian grid
TEST_CASE("PhaseProperties", "Water")
{
    auto water{FluidFactory::create_water(0.0, 0.0)};
}
