#pragma once

#include <vector>
#include <cassert>

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

    } // Grids

} // GPN