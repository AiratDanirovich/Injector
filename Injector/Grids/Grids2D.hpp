#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <numbers>

#include <Eigen/Core>

#include <Injector/Declarations.h>
#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateSystem.hpp>
#include <Injector/Grids/ConcreteGrids.hpp>
#include <Injector/Grids/GridsFactory.hpp>

namespace GPN
{
    namespace Grids
    {
        /// @brief Two-dimensional grid
        /// @tparam Axes1 Type for grid in first [r] direction
        /// @tparam Axes2 Type for grid in second [z] direction
        template <typename CoordinateSystem_t>
            requires IStructuredGrid2D<CoordinateSystem_t>
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
            using FaceAreaAxes2 = Eigen::Array<RealType, 1, -1>;
            struct Point
            {
                RealType x, y;
            };

            const AxesGrid<Axes1> first_coord;
            const AxesGrid<Axes2> second_coord;

            template <typename Axes_t>
            const auto &coordinate() const
            {
                if constexpr (std::is_same<Axes_t, Axes1>::value)
                    return first_coord;

                if constexpr (std::is_same<Axes_t, Axes2>::value)
                    return second_coord;
            }

            StructuredGrid2D(
                const AxesGrid<Axes1> &first_coord,
                const AxesGrid<Axes2> &second_coord)
                : first_coord{first_coord},
                  second_coord{second_coord},
                  its_volumes(
                      first_coord.volumes().matrix() *
                      second_coord.volumes().transpose().matrix())
            {
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
                return its_volumes(id1, id2);
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
            using StructuredGrid2D<CylinderCoordinates>::coordinate;

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
            FaceAreaAxes2 set_axes2_area() const
            {
                return FaceAreaAxes1{first_coord.volumes() * TwoPI()};
            }
        };

        struct StructuredXYGrid2D : public StructuredGrid2D<Cartesian2DCoordinates>
        {
            using StructuredGrid2D<Cartesian2DCoordinates>::StructuredGrid2D;
        };

        struct CylinderGridFactory
        {
            static auto create(const Box &box, ptrdiff_t n1, ptrdiff_t n2)
            {
                return create(GPN::Grids::Factory::generate_dual_grid_stencils_uniform(box.axes1, n1),
                              GPN::Grids::Factory::generate_dual_grid_stencils_uniform(box.axes2, n2));
            }

            static auto create(const auto &z_stencils, const auto &r_stencils)
            {
                return create_cylinder_grid_2D_ptr(
                    z_stencils, r_stencils);
            }
            
            template<typename Refiner_t>
            static auto create(Refiner_t&& refiner, const auto &z_stencils, const auto &r_stencils)
            {
                return create_cylinder_grid_2D_ptr(
                    refiner,
                    z_stencils, r_stencils);
            }

        protected:
            static auto create_cartesian_grid_2D_ptr(ptrdiff_t n)
            {
                auto stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, n)};

                auto nodes{GridDual{stencils, CoordinateTypes::X{}}};

                auto x_grid{
                    AxesGrid<CoordinateTypes::X>{nodes}};

                auto y_grid{
                    AxesGrid<CoordinateTypes::Y>{nodes}};

                return std::make_shared<StructuredXYGrid2D>(x_grid, y_grid);
            }

            static auto create_cylinder_grid_2D_ptr(const Box &box, ptrdiff_t n1, ptrdiff_t n2)
            {
                auto z_stencils{Factory::generate_dual_grid_stencils_uniform(box.axes1, n1)};
                auto r_stencils{Factory::generate_dual_grid_stencils_uniform(box.axes2, n2)};
                return create_cylinder_grid_2D_ptr(z_stencils, r_stencils);
            }
            
            static auto create_cylinder_grid_2D_ptr(const auto &z_stencils, const auto &r_stencils)
            {
                auto z_grid{Factory::create_axes<CoordinateTypes::Z>(z_stencils)};
                auto r_grid{Factory::create_axes<CoordinateTypes::R_CylCoord>(r_stencils)};

                return std::make_shared<StructuredCylinderGrid2DAxisymmetric>(z_grid, r_grid);
            }
            
            template<typename Refiner_t>
            static auto create_cylinder_grid_2D_ptr(Refiner_t&& refiner, const auto &z_stencils, const auto &r_stencils)
            {
                auto z_grid{Factory::create_axes<CoordinateTypes::Z>(refiner, z_stencils)};
                auto r_grid{Factory::create_axes<CoordinateTypes::R_CylCoord>(r_stencils)};

                return std::make_shared<StructuredCylinderGrid2DAxisymmetric>(z_grid, r_grid);
            }
        };

    } // Grids
} // GPN