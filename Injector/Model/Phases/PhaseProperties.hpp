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
    struct JouleThomson : public SomeProperty
    {
        // change units 
        JouleThomson(RealType value)
        : SomeProperty{value*1e-5}
        {}
    };

    struct StationaryPhaseProperties
    {
        StationaryPhaseProperties(const StationaryPhaseProperties&) = default;
        StationaryPhaseProperties(
            Density density,
            SpecificHeatCapacity specific_heat_capacity,
            HeatConductivity heat_conductivity) noexcept
            : density{density},
              specific_heat_capacity{specific_heat_capacity},
              heat_conductivity{heat_conductivity},
              volumetric_heat_capacity{
                  specific_heat_capacity *
                  density}
        {
        }
        const RealType density;
        const RealType specific_heat_capacity;
        const RealType heat_conductivity;
        const RealType volumetric_heat_capacity;
    };

    struct PhaseProperties : public StationaryPhaseProperties
    {
        PhaseProperties(
            Viscosity viscosity,
            Density density,
            SpecificHeatCapacity specific_heat_capacity,
            HeatConductivity heat_conductivity) noexcept
            : StationaryPhaseProperties{
                  density, specific_heat_capacity,
                  heat_conductivity},
              viscosity{viscosity}
        {
        }
        const RealType viscosity;
    };

    struct PhasePropertiesJT
    : public PhaseProperties
    {
        PhasePropertiesJT( 
            const PhaseProperties& props,
            JouleThomson JT) noexcept
            : PhaseProperties{props},
              JT{JT}
        {
        }

        const RealType JT;
    };

    struct Water : public PhaseProperties
    {
        Water(const PhaseProperties &props) noexcept
            : PhaseProperties{props}
        {
        }
    };

    struct WaterJT : public PhasePropertiesJT
    {
        WaterJT(const PhasePropertiesJT &props) noexcept
            : PhasePropertiesJT{props}
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