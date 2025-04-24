#pragma once

#include <Injector/Declarations.h>

namespace GPN
{
    template<typename Axes1_t, typename Axes2_t>
    struct CoordinateSystem2D
    {
        static constexpr size_t Dim()
        {
            return 2ull;
        }
        using Axes1 = Axes1_t;
        using Axes2 = Axes2_t;
    };
    
    template<typename Axes1_t, typename Axes2_t, typename Axes3_t>
    struct CoordinateSystem3D
    {
        static constexpr size_t Dim()
        {
            return 3ull;
        }
        using Axes1 = Axes1_t;
        using Axes2 = Axes2_t;
        using Axes3 = Axes3_t;
    };
} // GPN
