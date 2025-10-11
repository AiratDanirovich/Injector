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
            using cMarginMap1D = Eigen::Map<const DualNodesContainer>;

        public:
            GridDualStencilsMap(
                const GridDualStencils &
                    grid_dual_stencils)
                : dual_nodes{grid_dual_stencils.dual_nodes}
            // : dual_nodes{
            //       grid_dual_stencils.dual_nodes.data() + start_margin,
            //       grid_dual_stencils.dual_nodes - start_margin}
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

            const DualNodesContainer dual_nodes;
        };

        template <typename CoordinateType_t, ptrdiff_t start_margin>
        struct AxesGridMap
        {
            using cMarginMap1D = Eigen::Map<const ControlVolumesContainer>;

            AxesGridMap(const AxesGrid<CoordinateType_t> &axes_grid)
                : axes_grid{axes_grid},
                  dual_stencils{axes_grid.dual_stencils},
                  control_volumes{
                      axes_grid.volumes().data() + start_margin,
                      axes_grid.volumes().size() - start_margin}
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

            const GridDualStencilsMap<start_margin> dual_stencils;
            const cMarginMap1D control_volumes;

        private:
            const AxesGrid<CoordinateType_t> &axes_grid;
        };

    } // Grids

} // GPN