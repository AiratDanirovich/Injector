#pragma once

#include <vector>
#include <assert>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    namespace Grids
    {
        struct RefinerVerticle
        {
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
                const auto &nodes{dual_nodes_stencils.dual_nodes};

                RealType top{nodes.head(1ll)(0ll)},
                    bot{nodes.tail(1ll)(0ll)};
                ptrdiff_t layers{std::ceil((bot - top) / step)};

                std::vector<RealType> buf;
                buf.reserve(layers + 1ll);

                // push the top node
                buf.push_back(top);
                for (auto i{1ll}; i < is_permeable.size(); ++i)
                {
                    if (is_permeable(i) == 1.0)
                        // perforated layer -- do nothing
                        buf.push_back(nodes(i));
                    else if (is_permeable(i) == 0.0)
                    {
                        // rocks -- refine grid
                        RealType l_top{nodes(i - 1ll)}, l_bot{nodes(i)};
                        RealType thickness{l_bot - l_top};
                        ptrdiff_t segm_nmbr{std::ceil(thickness / step)};
                        RealType local_step{thickness / segm_nmbr};

                        for (auto j{0ll}; j < segm_nmbr - 1ll; ++j)
                            buf.push_back(buf.back() + local_step);
                        bur.push_back(l_bot);
                    }
                    else
                    // is_permeable should take only {0.0, 1.0} values
                        assert(false);
                }

                // nodes are in increasing order
                assert(buf[0ull] == nodes(0ll));
                for (auto i{1ull}; i < buf.size(); ++i)
                    assert(buf[i - 1ull] < buf[i]);
            }

        protected:
            const LogValuesContainer &is_permeable;
            const RealType step;
        };

    } // Grids

} // GPN