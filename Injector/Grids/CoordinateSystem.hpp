#pragma once

#include <Injector/Declarations.h>

namespace GPN
{
    /// @brief Definition of types of coordinates in 2D system
    /// @tparam Axes1_t First coordinate (X,Y,Z,R etc)
    /// @tparam Axes2_t Second coordinate (X,Y,Z,R etc)
    template <typename Axes1_t, typename Axes2_t>
    struct CoordinateSystem2D
    {
        static constexpr size_t Dim()
        {
            return 2ull;
        }
        using Axes1 = Axes1_t;
        using Axes2 = Axes2_t;
    };

    /// @brief 2D Cartesian coordinates, X-Y plane
    struct Cartesian2DCoordinates : public
        CoordinateSystem2D<
            CoordinateTypes::X,
            CoordinateTypes::Y>{};

    /// @brief 2D cylinder coordinates, Z-R plane
    struct CylinderCoordinates : public
        CoordinateSystem2D<
            CoordinateTypes::Z,
            CoordinateTypes::R_CylCoord>{};

    /// @brief Definition of types of coordinates in 3D system
    /// @tparam Axes1_t First coordinate (X,Y,Z,R etc)
    /// @tparam Axes2_t Second coordinate (X,Y,Z,R etc)
    /// @tparam Axes3_t Third coordinate (X,Y,Z,R etc)
    template <typename Axes1_t, typename Axes2_t, typename Axes3_t>
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

    /// @brief 3D Cartesian coordinates, X-Y-Z space
    struct Cartesian3DCoordinates : public
        CoordinateSystem3D<
            CoordinateTypes::X,
            CoordinateTypes::Y,
            CoordinateTypes::Z>{};
} // GPN
