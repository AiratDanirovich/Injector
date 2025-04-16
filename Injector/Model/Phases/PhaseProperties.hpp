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
    { }; 
    struct Density : public SomeProperty
    { }; 
    struct SpecificHeatCapacity : public SomeProperty
    { }; 

    struct PhaseProperties
    {
        PhaseProperties(
            Viscosity viscosity,
            Density density,
            SpecificHeatCapacity mass_heat_capacity) noexcept
            : viscosity{viscosity}
            , density{density}
            , mass_heat_capacity{mass_heat_capacity}
        {
            volumetric_heat_capacity = 
                    mass_heat_capacity*
                    density;
        }
        Viscosity viscosity;
        Density density;
        SpecificHeatCapacity mass_heat_capacity;
        RealType volumetric_heat_capacity;
    //    RealType temperature;
    //    RealType pressure;
    };
    
    struct Water : public PhaseProperties
    {
        Water(const PhaseProperties& props) noexcept
            : PhaseProperties{props}
        {}
    };
    
    struct Oil : public PhaseProperties
    {
        Oil(const PhaseProperties& props) noexcept
            : PhaseProperties{props}
        {}
    };
    
    struct Gas : public PhaseProperties
    {
        Gas(const PhaseProperties& props) noexcept
            : PhaseProperties{props}
        {}
    };

}// GPN