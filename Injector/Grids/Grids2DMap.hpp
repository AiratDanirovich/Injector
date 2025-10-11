#pragma once

#include <Eigen/Core>

#include <Injector/Grids/Grids2D.hpp>

namespace GPN
{
    namespace Grids
    {

        template <ptrdiff_t left_margin, typename Grid2D_t>
        struct Grid2DMap
        {
            Grid2DMap(const cptr<Grid2D_t> grid2D)
            :grid2D{grid2D}
            {}

        private:
            const cptr<Grid2D_t> grid2D;
        };

        using CylinderGridRock =
            Grid2DMap<3ll, StructuredCylinderGrid2DAxisymmetric>;

    } // Grids

} // GPN