#pragma once

#include <Eigen/Core>

#include <Injector/Grids/GridMap/Grids1DMap.hpp>
#include <Injector/Grids/Grids2D.hpp>

namespace GPN
{
    namespace Grids
    {
        template<ptrdiff_t left_margin, typename Grid2D_t>
        struct Grid2DMap
        {
            Grid2DMap(const cptr<Grid2D_t> grid2D)
            :grid2D{grid2D}
            {}

        private:
            using Axes1 = typename Grid2D_t::Axes_1;
            using Axes2 = typename Grid2D_t::Axes_2;

            const AxesGrid<Axes1> its_first_coord;
            const AxesGridMap<left_margin, Axes2> its_second_coord;
            const cptr<Grid2D_t> grid2D;
        };

        using CylinderGridRock =
            Grid2DMap<3ll, StructuredCylinderGrid2DAxisymmetric>;

    } // Grids

} // GPN