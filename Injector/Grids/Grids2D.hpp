#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

#include <Eigen/Core>
// #include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids1D.hpp>
#include <Injector/Grids/CoordinateTypes.h>

namespace GPN
{
    namespace Grids
    {
        /// @brief Two-dimensional grid
        /// @tparam FirstDir Type for grid in first [r] direction
        /// @tparam SecondDir Type for grid in second [z] direction
        template <typename FirstDir, typename SecondDir>
        struct StructuredGrid2D
        {
        public:
            struct Point
            {
                RealType x, y;
            };

            FirstDir first_coord;
            SecondDir second_coord;

            StructuredGrid2D(
                const FirstDir &first_coord, 
                const SecondDir &second_coord)
                    : first_coord{first_coord}
                    , second_coord{second_coord}
                    , its_volumes(
                        first_coord.size(),
                        second_coord.size()
                    )
            {
                // set volumes
                for (std::ptrdiff_t j = 0; j < its_volumes.cols(); ++j)
                    for (std::ptrdiff_t i = 0; i < its_volumes.rows(); ++i)
                        its_volumes(i, j) =
                            first_coord.volume(i) *
                            second_coord.volume(j);
            }

            auto ccordinates(auto id1, auto id2) const
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

        struct StructuredCylinderGrid2D
            : public StructuredGrid2D<
                  Grid1D<CoordinateTypes::Z>,
                  Grid1D<CoordinateTypes::RadialCylinderCoordinate>>
        {
            using StructuredGrid2D<
                Grid1D<CoordinateTypes::Z>,
                Grid1D<CoordinateTypes::RadialCylinderCoordinate>>::StructuredGrid2D;
        };

    } // Grids
} // GPN