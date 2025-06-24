#pragma once

#include <Injector/Grids/Defines.h>

namespace GPN
{
    struct SomeProperty
    {
        operator RealType() const { return value; }
        RealType value;
    };
    struct Viscosity : public SomeProperty
    {
    };
    struct Density : public SomeProperty
    {
    };
    struct SpecificHeatCapacity : public SomeProperty
    {
    };
    struct HeatConductivity : public SomeProperty
    {
    };

    struct StationaryPhaseProperties
    {
        StationaryPhaseProperties(const StationaryPhaseProperties&) = default;
        StationaryPhaseProperties(
            Density density,
            SpecificHeatCapacity mass_heat_capacity,
            HeatConductivity heat_conductivity) noexcept
            : density{density},
              mass_heat_capacity{mass_heat_capacity},
              heat_conductivity{heat_conductivity},
              volumetric_heat_capacity{
                  mass_heat_capacity *
                  density}
        {
        }
        const RealType density;
        const RealType mass_heat_capacity;
        const RealType heat_conductivity;
        const RealType volumetric_heat_capacity;
    };

    struct PhaseProperties : public StationaryPhaseProperties
    {
        PhaseProperties(
            Viscosity viscosity,
            Density density,
            SpecificHeatCapacity mass_heat_capacity,
            HeatConductivity heat_conductivity) noexcept
            : StationaryPhaseProperties{
                  density, mass_heat_capacity,
                  heat_conductivity},
              viscosity{viscosity}
        {
        }
        const RealType viscosity;
        //    RealType temperature;
        //    RealType pressure;
    };

    struct Water : public PhaseProperties
    {
        Water(const PhaseProperties &props) noexcept
            : PhaseProperties{props}
        {
        }
    };

    struct Oil : public PhaseProperties
    {
        Oil(const PhaseProperties &props) noexcept
            : PhaseProperties{props}
        {
        }
    };

    struct Gas : public PhaseProperties
    {
        Gas(const PhaseProperties &props) noexcept
            : PhaseProperties{props}
        {
        }
    };

} // GPN