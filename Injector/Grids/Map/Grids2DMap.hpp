#pragma once

#include <Eigen/Core>

#include <Injector/Grids/CoordinateSystem.hpp>
#include <Injector/Grids/Map/Grids1DMap.hpp>
#include <Injector/Grids/Grids2D.hpp>

namespace GPN
{
    namespace Grids
    {
        template <ptrdiff_t left_margin, typename Grid2D_t>
        struct Grid2DMap : public CoordinateSystem2D<typename Grid2D_t::Axes1, typename Grid2D_t::Axes2>
        {
            using Base = CoordinateSystem2D<typename Grid2D_t::Axes1, typename Grid2D_t::Axes2>;
            using typename Base::Axes1;
            using typename Base::Axes2;

            Grid2DMap(const cptr<Grid2D_t> grid2D)
                : grid2D{grid2D},
                its_volumes{grid2D->volumes().data(), 10ll, 10ll},
                its_first_coord{grid2D->first_coord()},
                its_second_coord{AxesGridMap<Axes2, left_margin>{grid2D->second_coord()}}
            {
            }

            auto volume(const auto id1, const auto id2) const
            {
                return volumes()(id1, id2);
            }
            const auto &volumes() const
            {
                return its_volumes;
            }

            const auto &first_coord() const
            {
                return its_first_coord;
            }
            const auto &second_coord() const
            {
                return its_second_coord;
            }

        private:
            using cMarginMap = Eigen::Map<const CellVolumeContainer2D>;

            const AxesGrid<Axes1>& its_first_coord;
            const AxesGridMap<Axes2, left_margin> its_second_coord;
            const cptr<Grid2D_t> grid2D;
            const cMarginMap its_volumes;
        };

        using CylinderGridRock =
            Grid2DMap<3ll, StructuredCylinderGrid2DAxisymmetric>;

    } // Grids

} // GPN