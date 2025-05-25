#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Phases;

RealType viscosity{6e-4}, density{1000}, capacity{4200};

TEST_CASE("PhaseProperties", "Water")
{
    auto water1{FluidFactory::create_water(0.0, 0.0)};
    const auto water2{
        FluidFactory::create_water(
            Viscosity{viscosity},
            Density{density},
            SpecificHeatCapacity{capacity})};
}
