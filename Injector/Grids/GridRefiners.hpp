#pragma once

#include <vector>
#include <cassert>
#include <algorithm>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Logs.hpp>

#include <Injector/Model/Well/WellHoles.hpp>

namespace GPN
{
    namespace Grids
    {
        struct EmptyRefiner
        {
            DualNodesContainer refine(
                const GridDualStencils &dual_nodes_stencils) noexcept
            {
                return refine(dual_nodes_stencils.dual_nodes);
            }
            DualNodesContainer refine(
                const DualNodesContainer &dual_nodes) noexcept
            {
                return dual_nodes;
            }

            DualNodesContainer refine(
                const std::vector<RealType> &buf) noexcept
            {
                DualNodesContainer out(buf.size());
                std::copy(buf.cbegin(), buf.cend(), out.begin());
                return out;
            }
        };

        struct RefinerVerticle
        {
            RefinerVerticle(
                const RealType step,
                LogValuesContainer &&is_permeable)
                : is_permeable{std::move(is_permeable)},
                  step{step}
            {
            }
            RefinerVerticle(
                const RealType step,
                const LogValuesContainer &is_permeable)
                : is_permeable{is_permeable},
                  step{step}
            {
            }

            DualNodesContainer refine(
                const GridDualStencils &dual_nodes_stencils) noexcept
            {
                return refine(dual_nodes_stencils.dual_nodes);
            }
            DualNodesContainer refine(
                const DualNodesContainer &dual_nodes) noexcept
            {
                std::vector<RealType> buf(dual_nodes.size());
                std::copy(dual_nodes.cbegin(), dual_nodes.cend(), buf.begin());

                return refine(buf);
            }

            DualNodesContainer refine(
                const std::vector<RealType> &nodes) noexcept
            {
                for (auto i{0ll}; i < is_permeable.size(); ++i)
                    assert((is_permeable(i) == 0.0) || (is_permeable(i) == 1.0));

                assert(nodes.size() == is_permeable.size() + 1ll);

                RealType top{nodes.front()},
                    bot{nodes.back()};
                ptrdiff_t layers{static_cast<ptrdiff_t>(std::ceil((bot - top) / step))};

                std::vector<RealType> buf;
                buf.reserve(layers + 1ll);

                // push the top node
                buf.push_back(top);
                for (auto i{0ll}; i < is_permeable.size(); ++i)
                {
                    if (is_permeable(i) == 1.0)
                        // perforated layer -- do nothing
                        buf.push_back(nodes[i + 1ull]);
                    else if (is_permeable(i) == 0.0)
                    {
                        // rocks -- refine grid
                        RealType l_top{nodes[i]}, l_bot{nodes[i + 1ll]};
                        RealType thickness{l_bot - l_top};
                        ptrdiff_t segm_nmbr{static_cast<ptrdiff_t>(std::ceil(thickness / step))};
                        RealType local_step{thickness / segm_nmbr};

                        for (auto j{0ll}; j < segm_nmbr - 1ll; ++j)
                            buf.push_back(buf.back() + local_step);
                        buf.push_back(l_bot);
                    }
                    else
                        // is_permeable should take only {0.0, 1.0} values
                        assert(false);
                }

                // nodes are in increasing order
                assert(buf[0ull] == nodes[0ull]);
                for (auto i{1ull}; i < buf.size(); ++i)
                    assert(buf[i - 1ull] < buf[i]);

                DualNodesContainer out(buf.size());
                std::copy(buf.cbegin(), buf.cend(), out.begin());
                return out;
            }

        protected:
            const LogValuesContainer is_permeable;
            const RealType step;
        };

        struct RefinerRadial_LogWellHoles
        {
            RefinerRadial_LogWellHoles(
                const WellHoles &well_holes,
                const RealType rMax,
                const RealType q,
                const RealType max_step)
                : tube_inner_radius{well_holes.tube_inner_radius},
                  sandface_radius{well_holes.sandface_radius},
                  column_outer_radius{well_holes.column_outer_radius},
                  r_max{rMax},
                  q{q},
                  max_step{max_step},
                  r_min{0.0}
            {
            }

            DualNodesContainer refine(
                const GridDualStencils &dual_nodes_stencils) noexcept
            {
                return refine(dual_nodes_stencils.dual_nodes);
            }
            DualNodesContainer refine(
                const DualNodesContainer &dual_nodes) noexcept
            {
                std::vector<RealType> buf(dual_nodes.size());
                std::copy(dual_nodes.cbegin(), dual_nodes.cend(), buf.begin());

                return refine(buf);
            }

            DualNodesContainer refine(
                const std::vector<RealType> &nodes) noexcept
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

                    DualNodesContainer buf(out.size());
                    std::copy(out.cbegin(), out.cend(), buf.begin());
                    return buf;
                }
            }

        protected:
            const RealType sandface_radius;
            const RealType tube_inner_radius;
            const RealType column_outer_radius;
            const RealType r_max;
            const RealType q;
            const RealType max_step;
            const RealType r_min;

        private:
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

            DualNodesContainer generate_uniform_radial_grid(
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

                DualNodesContainer buf(out.size());
                std::copy(out.cbegin(), out.cend(), buf.begin());

                return buf;
            }
        };
    } // Grids
} // GPN