#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>

#include <Eigen/Core>
// #include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateSystem.hpp>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids1D.hpp>

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
            //    static_assert(std::is_same<CoordinateSystem_t::Axes1, CoordinateTypes::Z>::value);
            //    static_assert(std::is_same<CoordinateSystem_t::Axes2, CoordinateTypes::R_CylCoord>::value);

            using typename CoordinateSystem_t::Axes1;
            using typename CoordinateSystem_t::Axes2;

        public:
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

        using CylinderCoordinates =
            CoordinateSystem2D<
                CoordinateTypes::Z,
                CoordinateTypes::R_CylCoord>;
        using Cartesian2DCoordinates =
            CoordinateSystem2D<
                CoordinateTypes::X,
                CoordinateTypes::Y>;

        using StructuredCylinderGrid2D =
            StructuredGrid2D<CylinderCoordinates>;

    } // Grids
} // GPN