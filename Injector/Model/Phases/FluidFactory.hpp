#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>

namespace GPN
{
    namespace Phases
    {
        struct FluidFactory
        {
            static auto create_water(RealType temperature, RealType pressure)
            {
                return Water{
                    PhaseProperties{
                        Viscosity{6e-4},
                        Density{1000},
                        SpecificHeatCapacity{4180}}};
            }
            static auto create_water(Viscosity viscosity,
                                     Density density,
                                     SpecificHeatCapacity mass_heat_capacity)
            {
                return Water{
                    PhaseProperties{
                        viscosity,
                        density,
                        mass_heat_capacity}};
            }
        };

    } // Phases
} // GPN
