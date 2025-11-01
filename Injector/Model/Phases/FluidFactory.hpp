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
                         SpecificHeatCapacity mass_heat_capacity,
                         HeatConductivity heat_conductivity)
                : water{create_water(
                      viscosity,
                      density,
                      mass_heat_capacity,
                      heat_conductivity)}
            {
            }

            PhaseProperties water;
            static Water create_water(RealType temperature, RealType pressure)
            {
                return {
                    PhaseProperties{
                        Viscosity{6e-4},
                        Density{1000},
                        SpecificHeatCapacity{4000},
                        HeatConductivity{0.6}}};
            }

            static Water create_water(
                Viscosity viscosity,
                Density density,
                SpecificHeatCapacity mass_heat_capacity,
                HeatConductivity heat_conductivity)
            {
                return {
                    PhaseProperties{
                        viscosity,
                        density,
                        mass_heat_capacity,
                        heat_conductivity}};
            }

            static WaterJT create_water_JT(
                Viscosity viscosity,
                Density density,
                SpecificHeatCapacity mass_heat_capacity,
                HeatConductivity heat_conductivity,
                JouleThomson joule_thomson,
                const RealType adiabatic_weight)
            {
                return WaterJT{
                    PhasePropertiesJT{
                        PhaseProperties{
                            viscosity,
                            density,
                            mass_heat_capacity,
                            heat_conductivity},
                        joule_thomson,
                        adiabatic_weight}};
            }
        };

    } // Phases
} // GPN
