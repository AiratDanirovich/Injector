#pragma once

#include <cmath>

#include "defines.h"

namespace GPN
{
    namespace CoordinateTypes
    {

        struct CartesianCoordinate
        {
            static float_t volume(float_t x1, float_t x2)
            {
                return x2 - x1;
            }
            static float_t resistivity(float_t x1, float_t x2)
            {
                return x2 - x1;
            }
        };

        struct RadialCylinderCoordinate
        {
            static float_t volume(float_t x1, float_t x2)
            {
                assert(x2 > x1);
                return (x2 * x2 - x1 * x1) / 2.0;
            }
            static float_t resistivity(float_t x1, float_t x2)
            {
                assert(x2 > x1);
                assert(x1 > 0.0);
                return std::log(x2/x1);
            }
        };
    } // CoordinateTypes
} // GPN