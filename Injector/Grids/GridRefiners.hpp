#pragma once

#include <vector>
#include <cassert>
#include <algorithm>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Logs.hpp>

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

        struct RefinerVerticle_CrossFlow
        {
            RefinerVerticle_CrossFlow(
                const RealType regular_step,
                const LogValuesContainer &is_permeable,
                const std::vector<RealType> &from_coord,
                const RealType crossflow_step)
                : is_permeable{is_permeable},
                  regular_step{regular_step},
                  from_coord{from_coord},
                  crossflow_step{crossflow_step}
            {
                std::sort(this->from_coord.begin(), this->from_coord.end());

                const auto& this_coord{this->from_coord};
                for (auto i{1ull}; i < this_coord.size(); ++i)
                    assert(this_coord[i] >= this_coord[i - 1ull] + min_spasing);
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
                const std::vector<RealType> &dual_nodes) noexcept
            {
                std::vector<RealType> buf;
                buf.reserve(dual_nodes.size() + from_coord.size() + 2ull);

                RealType top{dual_nodes.front()},
                    bot{dual_nodes.back()};
                ptrdiff_t layers{static_cast<ptrdiff_t>(std::ceil((bot - top) / regular_step))};

                auto cf_id{0ull};
                // push the top node
                buf.push_back(top);
                for (auto i{0ll}; i < is_permeable.size(); ++i)
                {
                    const RealType l_top{dual_nodes[i]};
                    const RealType l_bot{dual_nodes[i + 1ll]};
                    if (is_permeable(i) == 1.0)
                    {
                        // perforated layer -- do nothing
                        buf.push_back(l_bot);
                        // skip cross_flow coords above the current layer bottom
                        while (from_coord[cf_id] < l_bot)
                            ++cf_id;
                    }
                    else if (is_permeable(i) == 0.0)
                    {
                        // rocks -- refine grid
                        // treat the cross_flow coords which are within the current layer
                        while (
                            (cf_id < from_coord.size()) &&
                            (from_coord[cf_id] < l_bot))
                        {
                            const RealType cf_top{std::max(l_top, from_coord[cf_id] - min_spacing / 2.0)};
                            assert(cf_top >= l_top);
                            assert(cf_top < l_bot);
                            assert(buf.back() <= cf_top);
                            if (buf.back() < cf_top)
                            {
                                // rock layer top is above the cross-flow sublayer top.
                                // refine rock sublayer and push the refined segments into the buf
                                insert_rock_sublayers(buf, buf.back(), cf_top);
                            }

                            const RealType cf_bot{cf_top + min_spacing};
                            assert(cf_bot <= l_bot);
                            assert(cf_bot > l_top);
                            assert(buf.back() < cf_bot);
                            buf.push_back(cf_bot);

                            ++cf_id;
                        }
                        // push bottom boundary only
                        // if it is below the current cross_flow coordinate
                        if (buf.back() < l_bot)
                            insert_rock_sublayers(buf, buf.back(), l_bot);
                    }
                    else
                    {
                        // is_permeable should take only {0.0, 1.0} values
                        assert(false);
                    }
                }

                // nodes must be in STRICTLY increasing order
                for (auto i{1ull}; i < buf.size(); ++i)
                    assert(buf[i - 1ull] < buf[i]);
                assert(buf.back() == nodes.back());

                DualNodesContainer out(buf.size());
                std::copy(buf.cbegin(), buf.cend(), out.begin());
                return out;
            }

        protected:
            const LogValuesContainer is_permeable;
            const RealType regular_step;
            std::vector<RealType> from_coord;
            const RealType crossflow_step;

            // minimum distance between two
            // cross_flow coordinates
            const RealType min_spacing{1.0};

        private:
            void insert_rock_sublayers(std::vector<RealType> &buf, const RealType top, const RealType bot)
            {
                // insert rock sublayers between [top < bot]
                const RealType thickness{bot - top};
                assert(thickness > 0.0);
                // nmbr of segments to insert
                const ptrdiff_t segm_nmbr{static_cast<ptrdiff_t>(std::ceil(thickness / regular_step))};
                const RealType local_step{thickness / segm_nmbr};
                for (auto j{0ll}; j < segm_nmbr - 1ll; ++j)
                    buf.push_back(buf.back() + local_step);
                // push cf top as bottom of Rock sublayer
                buf.push_back(bot);
            }
        };

        struct AbstractRefinerRadial
        {
            DualNodesContainer refine(
                const GridDualStencils &dual_nodes_stencils) const noexcept
            {
                return refine(dual_nodes_stencils.dual_nodes);
            }
            DualNodesContainer refine(
                const DualNodesContainer &dual_nodes) const noexcept
            {
                std::vector<RealType> buf(dual_nodes.size());
                std::copy(dual_nodes.cbegin(), dual_nodes.cend(), buf.begin());

                return refine(buf);
            }

            virtual DualNodesContainer refine(
                const std::vector<RealType> &r_stencils) const noexcept = 0;

        protected:
            std::vector<RealType> init_grid(
                const std::vector<RealType> &r_stencils,
                const std::ptrdiff_t r_nodes) const
            {
                assert(r_stencils.size() >= 2ll);
                std::vector<RealType> out;
                out.reserve(r_nodes + 5ull);
                // forget the r_max
                for (auto i{0ull}; i < r_stencils.size() - 1ull; ++i)
                    out.push_back(r_stencils[i]);
                return out;
            }

            DualNodesContainer generate_uniform_radial_grid(
                const std::vector<RealType> &r_stencils,
                const ptrdiff_t r_nodes) const
            {
                assert(r_stencils.size() >= 2ull);

                const RealType r_max{r_stencils.back()};
                const RealType r_min{r_stencils.front()};
                const auto sandface_radius{r_stencils[r_stencils.size() - 2ll]};

                assert(r_min == 0.0);
                assert(r_max > sandface_radius);
                assert(r_nodes > 1ll);

                // forget r_max in the origianl stencils container
                std::vector<RealType> out{init_grid(r_stencils, r_nodes)};
                // step of uniform mesh
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

        struct RefinerRadial_UniformWellHoles
            : public AbstractRefinerRadial
        {
            RefinerRadial_UniformWellHoles(
                const std::ptrdiff_t r_nodes_nmbr)
                : r_nodes_nmbr{r_nodes_nmbr}
            {
            }

            /// @brief
            /// @param r_stencils Includes r = 0, radii of sandwich materials, and r_max
            /// @return
            DualNodesContainer refine(
                const std::vector<RealType> &r_stencils) const noexcept override
            {
                assert(r_stencils.size() == 5ull);

                const auto r_min{r_stencils.front()},
                    tube_inner_radius{r_stencils[1ull]},
                    column_outer_radius{r_stencils[2ull]},
                    sandface_radius{r_stencils[3ull]},
                    r_max{r_stencils.back()};

                assert(r_min == 0.0);
                assert(r_min < tube_inner_radius);
                assert(tube_inner_radius < column_outer_radius);
                assert(column_outer_radius < sandface_radius);
                assert(sandface_radius < r_max);

                // uniform grid
                return generate_uniform_radial_grid(
                    r_stencils, r_nodes_nmbr);
            }

        protected:
            const std::ptrdiff_t r_nodes_nmbr;
        };

        struct RefinerRadial_LogWellHoles
            : public AbstractRefinerRadial
        {
            RefinerRadial_LogWellHoles(
                const RealType q,
                const RealType max_step)
                : q{q},
                  max_step{max_step}
            {
                assert(q >= 1.0);
            }

            /// @brief
            /// @param r_stencils Includes r = 0, radii of sandwich materials, and r_max
            /// @return
            DualNodesContainer refine(
                const std::vector<RealType> &r_stencils) const noexcept override
            {
                assert(r_stencils.size() == 5ull);

                const auto r_min{r_stencils.front()},
                    tube_inner_radius{r_stencils[1ull]},
                    column_outer_radius{r_stencils[2ull]},
                    sandface_radius{r_stencils[3ull]},
                    r_max{r_stencils.back()};

                assert(r_min == 0.0);
                assert(r_min < tube_inner_radius);
                assert(tube_inner_radius < column_outer_radius);
                assert(column_outer_radius < sandface_radius);
                assert(sandface_radius < r_max);
                // base, minimum step for geometric progression
                const auto base_step{sandface_radius - column_outer_radius};

                if (q == 1.0)
                {
                    // uniform grid
                    return generate_uniform_radial_grid(
                        r_stencils,
                        static_cast<ptrdiff_t>(
                            std::ceil((r_max - r_min) / base_step)));
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

                    std::vector<RealType> out{init_grid(r_stencils, nx)};

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
                    assert(out.back() >= r_max);

                    DualNodesContainer buf(out.size());
                    std::copy(out.cbegin(), out.cend(), buf.begin());
                    return buf;
                }
            }

        protected:
            const RealType q;
            const RealType max_step;
        };
    } // Grids
} // GPN