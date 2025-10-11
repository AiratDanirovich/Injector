#pragma once

#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace Grids
    {
        template <typename CoordinateType_t, ptrdiff_t start_margin>
        struct AxesGridMap
        {
            AxesGridMap(const AxesGrid<CoordinateType_t> &axes_grid)
            : axes_grid{axes_grid}
            {
            }

            auto mesh_size() const
            {
                return axes_grid.mesh_size() - start_margin;
            }

        private:
            const AxesGrid<CoordinateType_t> &axes_grid;
        };

    } // Grids

} // GPN