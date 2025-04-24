#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <numbers>

#include <Eigen/Core>
// #include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateSystem.hpp>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/ConcreteGrids.hpp>

namespace GPN
{
    namespace Grids
    {
        /// @brief Two-dimensional grid
        /// @tparam Axes1 Type for grid in first [r] direction
        /// @tparam Axes2 Type for grid in second [z] direction

        template <typename CoordinateSystem_t>
        struct StructuredGrid2D : public CoordinateSystem_t
        {
            static_assert(
                !std::is_same<
                    typename CoordinateSystem_t::Axes1,
                    typename CoordinateSystem_t::Axes2>::value);

            using typename CoordinateSystem_t::Axes1;
            using typename CoordinateSystem_t::Axes2;
            using CoordinateSystem = CoordinateSystem_t;

        public:
            using FaceAreaAxes1 = Eigen::ArrayX<RealType>;
            using FaceAreaAxes2 = Eigen::ArrayX<RealType>;
            struct Point
            {
                RealType x, y;
            };

            const AxesGrid<Axes1> first_coord;
            const AxesGrid<Axes2> second_coord;

            StructuredGrid2D(
                const AxesGrid<Axes1> &first_coord,
                const AxesGrid<Axes2> &second_coord)
                : first_coord{first_coord}, second_coord{second_coord}, its_volumes(
                                                                            first_coord.size(),
                                                                            second_coord.size())
            {
                // set volumes
                for (std::ptrdiff_t j = 0; j < its_volumes.cols(); ++j)
                    for (std::ptrdiff_t i = 0; i < its_volumes.rows(); ++i)
                        its_volumes(i, j) =
                            first_coord.volume(i) *
                            second_coord.volume(j);
            }

            auto coordinates(auto id1, auto id2) const
            {
                return Point{first_coord.mesh_nodes(id1), second_coord.mesh_nodes(id2)};
            }

            // steps in two directions,
            // between nodes id1 and id1+1 in first direction
            // and between nodes id2 and id2+1 in second direction
            auto steps(auto id1, auto id2) const
            {
                return {first_coord.mesh_steps(id1), second_coord.mesh_steps(id2)};
            }

            /// @brief cell volume at node ids {id1, id2}
            auto volume(auto id1, auto id2) const
            {
                return first_coord.cell_volumes(id1) * second_coord.cell_volumes(id2);
            }

            const auto &volumes() const
            {
                return its_volumes;
            }

            CellVolumeContainer2D its_volumes;
        };

        // take the third -- phi -- axes into account and multiply
        // all face areas by 2PI
        struct StructuredCylinderGrid2DAxisymmetric
            : public StructuredGrid2D<CylinderCoordinates>
        {
            constexpr static auto TwoPI()
            {
                return static_cast<RealType>(2.0 * std::numbers::pi);
            }
            StructuredCylinderGrid2DAxisymmetric(
                const AxesGrid<Axes1> &first_coord,
                const AxesGrid<Axes2> &second_coord)
                : StructuredGrid2D<CylinderCoordinates>{first_coord, second_coord},
                  face_area_axes1{set_axes1_area()},
                  face_area_axes2{set_axes2_area()}
            {
                // take axial symmetry into account,
                // multiply 2D-volumes by 2Pi
                its_volumes = its_volumes * TwoPI();
            }

            const FaceAreaAxes1 face_area_axes1;
            const FaceAreaAxes2 face_area_axes2;

        protected:
            FaceAreaAxes1 set_axes1_area() const
            {
                return FaceAreaAxes1{second_coord.volumes() * TwoPI()};
            }
            FaceAreaAxes1 set_axes2_area() const
            {
                return FaceAreaAxes1{first_coord.volumes() * TwoPI()};
            }
        };

        using StructuredXYGrid2D =
            StructuredGrid2D<Cartesian2DCoordinates>;

    } // Grids
} // GPN