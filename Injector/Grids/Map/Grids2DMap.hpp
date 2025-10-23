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

        private:
            using cMarginMap2D = Eigen::Map<const CellVolumeContainer2D>;
            using cMarginMap1D = Eigen::Map<const typename Grid2D_t::FaceAreaAxes1>;
            using cMarginMap1D_T = Eigen::Map<const typename Grid2D_t::FaceAreaAxes2>;

        public:
            using Base = CoordinateSystem2D<typename Grid2D_t::Axes1, typename Grid2D_t::Axes2>;
            using typename Base::Axes1;
            using typename Base::Axes2;

            const ptrdiff_t l_margin{left_margin};

            Grid2DMap(const cptr<Grid2D_t> grid2D)
                : grid2D{grid2D},
                  its_volumes{
                      grid2D->volumes().data() + left_margin * grid2D->first_coord().mesh_size(),
                      grid2D->first_coord().mesh_size(),
                      grid2D->second_coord().mesh_size() - left_margin},
                  its_first_coord{grid2D->first_coord()},
                  its_second_coord{AxesGridMap<Axes2, left_margin>{grid2D->second_coord()}},
                  face_area_axes1{grid2D->face_area_axes1.data() + left_margin, grid2D->face_area_axes1.size() - left_margin},
                  face_area_axes2{grid2D->face_area_axes2}
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

            template <typename Axes_t>
            const auto &coordinate() const
            {
                if constexpr (std::is_same<Axes_t, Axes1>::value)
                    return first_coord();

                if constexpr (std::is_same<Axes_t, Axes2>::value)
                    return second_coord();
            }
            
            /// @brief Rowmajor linear enumeration of mesh nodes
            /// @param first Node index along the first coordinate
            /// @param second Node index along the second coordinate
            /// @return Linear index, continuos numbering of all 2D nodes
            auto to_linear(const ptrdiff_t first, const ptrdiff_t second) const
            {
                assert(first >= 0ll);
                assert(second >= 0ll);
                assert(first < first_coord().mesh_size());
                assert(second < second_coord().mesh_size());
                return first + second * first_coord().mesh_size();
            }

            const auto &first_coord() const
            {
                return its_first_coord;
            }
            const auto &second_coord() const
            {
                return its_second_coord;
            }
            
            auto mesh_size() const
            {
                return first_coord().mesh_size() * second_coord().mesh_size();
            }

            using FaceAreaAxes1 = cMarginMap1D;
            using FaceAreaAxes2 = Grid2D_t::FaceAreaAxes2;
            const FaceAreaAxes1 face_area_axes1;
            const FaceAreaAxes2 &face_area_axes2;

        private:
            const AxesGrid<Axes1> &its_first_coord;
            const AxesGridMap<Axes2, left_margin> its_second_coord;
            const cptr<Grid2D_t> grid2D;
            const cMarginMap2D its_volumes;
        };

        using CylinderGridRock =
            Grid2DMap<3ll, StructuredCylinderGrid2DAxisymmetric>;

    } // Grids

} // GPN