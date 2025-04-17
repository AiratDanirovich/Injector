#pragma once

namespace GPN
{
    template<typename Axes1_t, typename Axes2_t>
    struct CoordinateSystem2D
    {
        using Axes1 = Axes1_t;
        using Axes2 = Axes2_t;
    };
    
    template<typename Axes1_t, typename Axes2_t, typename Axes2_t>
    struct CoordinateSystem3D : public CoordinateSystem2D<Axes1_t, Axes2_t>
    {
        using Axes3 = Axes3_t;
    };
} // GPN
