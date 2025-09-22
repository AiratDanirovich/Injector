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

        std::vector<RealType> generate_uniform_radial_grid(
            const RealType r_min,
            const RealType r_max,
            const ptrdiff_t r_nodes) const
        {
            assert(r_min < tube_inner_radius);
            assert(r_max > sandface_radius);
            assert(r_nodes > 1ll);

            std::vector<RealType> out{init_grid(r_min, r_max, r_nodes)};

            const RealType step{(r_max - sandface_radius) / (r_nodes - 1ll)};
            for (auto i{2ll}; i < r_nodes; ++i)
                out.push_back(out.back() + step);
            out.push_back(r_max);

            assert(out.front() == r_min);
            for (auto i{1ull}; i < out.size(); ++i)
                assert(out[i] > out[i - 1ull]);
            assert(out.back() == r_max);

            return out;
        }

        std::vector<RealType> generate_log_radial_grid(
            const RealType r_min, const RealType r_max,
            const RealType q,        // ratio of adjascent steps
            const RealType max_step) // max allowed step
            const
        {
            assert(r_min < tube_inner_radius);
            assert(r_max > sandface_radius);
            //    assert(max_step > tube_radius);
            //    assert(max_step > sandface_radius - tube_radius);
            assert(q >= 1.0);
            // base, minimum step for geometric progression
            const RealType base_step = sandface_radius - column_outer_radius;

            if (q == 1.0)
            {
                // uniform grid
                return generate_uniform_radial_grid(r_min, r_max, (ptrdiff_t)std::ceil((r_max - r_min) / base_step));
            }
            else
            {
                // non-uniform grid
                // nmbr of steps within segment [sandface_radius; r_max],
                // which increase geometrically
                ptrdiff_t nx{
                    (ptrdiff_t)std::ceil(
                        std::log(1.0 + (r_max - sandface_radius) / base_step * (q - 1.0)) /
                        std::log(q))};

                std::vector<RealType> out{init_grid(r_min, r_max, nx)};

                // recalculate the base step
                const RealType hx{base_step}; //{(r_max - sandface_radius) * (q - 1.0) / (std::pow(q, nx) - 1.0)};
                assert(hx <= base_step);
                assert(hx > 0.0);

                if (max_step < hx * (std::pow(q, nx - 1ll)))
                {
                    // the furthest steps are too large
                    for (auto i{0ll}; i < nx; ++i)
                        out.push_back(out.back() + std::min(max_step, hx * std::pow(q, i)));
                    // tail of the segment where the geometric steps are too large
                    while (out.back() < r_max)
                        out.push_back(out.back() + max_step);
                }
                else
                {
                    // the furthest steps are fine
                    for (auto i{0ll}; i < nx; ++i)
                        out.push_back(out.back() + hx * std::pow(q, i));
                }
                for (auto i{1ull}; i < out.size(); ++i)
                    assert(out[i] > out[i - 1ull]);

                return out;
            }
        }

        std::vector<RealType> init_grid(
            const RealType r_min,
            const RealType r_max,
            const ptrdiff_t r_nodes) const
        {
            std::vector<RealType> out;
            out.reserve(r_nodes + 3ull);

            out.push_back(r_min);               // push leftmost boundary
            out.push_back(tube_inner_radius);   // push tube radius
            out.push_back(column_outer_radius); // push tube radius
            out.push_back(sandface_radius);     // push sandface radius

            return out;
        }

        const RealType sandface_radius;
        const RealType tube_inner_radius;
        const RealType column_outer_radius;
    };
} // GPN