#pragma once

#include <vector>
#include <cassert>

#include <Injector/Grids/Defines.h>

namespace GPN
{
    struct TubeInnerRadius : public SomeProperty
    {
    };
    struct ColumnOuterRadius : public SomeProperty
    {
    };
    struct SandfaceRadius : public SomeProperty
    {
    };

    /// @brief Descriptor for the well circular desing,
    /// contains radii of tube < column < sandface.
    /// Cement is between column and sandface
    struct WellHoles
    {
        WellHoles(
            const TubeInnerRadius tube_inner_radius,
            const ColumnOuterRadius column_outer_radius,
            const SandfaceRadius sandface_radius)
            : tube_inner_radius{tube_inner_radius},
              sandface_radius{sandface_radius},
              column_outer_radius{column_outer_radius}
        {
            assert(tube_inner_radius < column_outer_radius);
            assert(column_outer_radius < sandface_radius);
        }

        auto get_stencils(
            const RealType r_min, const RealType r_max) const noexcept
        {
            return std::vector<RealType>{
                r_min, tube_inner_radius,
                column_outer_radius, sandface_radius,
                r_max};
        }

        const RealType sandface_radius;
        const RealType tube_inner_radius;
        const RealType column_outer_radius;
    };
} // GPN