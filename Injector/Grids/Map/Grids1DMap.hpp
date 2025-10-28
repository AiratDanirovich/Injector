#pragma once

#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace Grids
    {
        template <ptrdiff_t start_margin>
        struct GridDualStencilsMap
        {
        private:
            using cMarginMap1D = Eigen::Map<const DualNodesContainer::Base>;

        public:
            GridDualStencilsMap(
                const GridDualStencils &
                    grid_dual_stencils)
                : dual_nodes{
                      grid_dual_stencils.dual_nodes.data() + start_margin,
                      grid_dual_stencils.dual_nodes.size() - start_margin}
            {
            }

            const auto &get_dual_nodes() const
            {
                return dual_nodes;
            }

            auto operator()(auto idx) const
            {
                assert(static_cast<ptrdiff_t>(idx) < static_cast<ptrdiff_t>(dual_nodes.size()));
                return dual_nodes(idx);
            }

            auto size() const
            {
                return dual_nodes.size();
            }

            const cMarginMap1D dual_nodes;
        };

        template <typename CoordinateType_t, ptrdiff_t start_margin>
        struct AxesGridMap
        {
        private:
            using cMarginMap1D = Eigen::Map<const ControlVolumesContainer>;

        public:
            using Axes = CoordinateType_t;

            AxesGridMap(const AxesGrid<CoordinateType_t> &axes_grid)
                : axes_grid{axes_grid},
                  control_volumes{
                      axes_grid.volumes().data() + start_margin,
                      axes_grid.volumes().size() - start_margin},
                  mesh_nodes{
                      axes_grid.mesh_nodes.data() + start_margin,
                      axes_grid.mesh_nodes.size() - start_margin},
                  dual_nodes{
                      axes_grid.dual_nodes.data() + start_margin,
                      axes_grid.dual_nodes.size() - start_margin}
            {
            }

            auto mesh_size() const
            {
                return axes_grid.mesh_size() - start_margin;
            }
            auto dual_size() const
            {
                return axes_grid.dual_size() - start_margin;
            }

            auto dual_front() const
            {
                return dual_nodes(0ll);
            }
            auto dual_back() const
            {
                return dual_nodes(dual_nodes.size() - 1ll);
            }

            const cMarginMap1D control_volumes;

            const cMarginMap1D mesh_nodes;
            // Dual mesh to be used in simulation
            const cMarginMap1D dual_nodes;

        private:
            const AxesGrid<CoordinateType_t> &axes_grid;
        };

    } // Grids

} // GPN