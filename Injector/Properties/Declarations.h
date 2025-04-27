#pragma once

#include <Injector/Grids/Declarations.h>

namespace GPN
{
    namespace Logs
    {

        struct InterpolatedDataContainer;
        struct StepProperty;
        struct StepPropertyGrid;
        struct AssertNonNegative;

        struct IndicatorProperty;
        struct IsPermeable;

        struct Permeability;
        struct Porosity;
        struct SkinFactor;
        struct SolidDensity;
        struct SolidSpecificHeatCapacity;
        struct SolidVolumetricHeatCapacity;
        struct HeatVolumetricCapacity;

        template <typename CoordinateType_t /* = CoordinateTypes::Z*/>
        struct FaceInterpolator;

        template <typename CoordinateType_t /* = CoordinateTypes::Z*/>
        struct FaceInterpolatedProperty;

        using ZInterpolator =
            FaceInterpolatedProperty<CoordinateTypes::Z>;

            struct HeatConductivity;
            struct ThermalDiffusivity;


    } // Logs
} // GPN

