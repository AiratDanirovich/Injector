#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>

namespace GPN
{
    namespace Phases
    {
        struct FluidFactory
        {
            FluidFactory(RealType temperature, RealType pressure)
                : water{create_water(temperature, pressure)}
            {
            }
            FluidFactory(Viscosity viscosity,
                         Density density,
                         SpecificHeatCapacity mass_heat_capacity)
                : water{create_water(viscosity,
                                     density,
                                     mass_heat_capacity)}
            {
            }

            PhaseProperties water;
            static Water create_water(RealType temperature, RealType pressure)
            {
                return {
                    PhaseProperties{
                        Viscosity{6e-4},
                        Density{1000},
                        SpecificHeatCapacity{4180}}};
            }
            static Water create_water(Viscosity viscosity,
                                     Density density,
                                     SpecificHeatCapacity mass_heat_capacity)
            {
                return {
                    PhaseProperties{
                        viscosity,
                        density,
                        mass_heat_capacity}};
            }
        };

    } // Phases
} // GPN
