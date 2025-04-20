#pragma once

#include <type_traits>
#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    namespace Logs
    {
        template <>
        struct FaceInterpolator<CoordinateTypes::R_CylCoord>
        {
            using Axes = CoordinateTypes::R_CylCoord;

            // const Properties::FaceAxes2VauesContainer face_values;

            // template <typename Grid_t>
            // FaceInterpolator(
            //     const StepPropertyGrid &log,
            //     const Grid_t &grid)
            //     : face_values{
            //           interpolate(log, grid)}
            // {
            //     static_assert(std::is_same<Axes, typename Grid_t::Axes>::value);
            // }

        public:
            template <typename Grid_t>
            static auto interpolate(
                const StepPropertyGrid &log,
                const Grid_t &grid)
            {
                // if dual_size() == 2 --- no internal faces
                assert(grid.dual_size() >= 2ll);

                FaceValuesContainer out(log.grid.mesh_size(), grid.dual_size() - 2ll);

                for (auto id{0ll}; id < out.cols(); ++id)
                    out.col(id) = Axes::const_face_interpolator(
                        grid.mesh_nodes(id), grid.mesh_nodes(id + 1ll),
                        static_cast<Logs::StepPropertyContainer>(log.log_vals));
                return out;
            }
        };
    } // Logs

    namespace Properties
    {
        /// @brief Interpolate physical property between nodes of 2D grid
        /// @tparam Grid_t 2D grid
        template <typename Grid_t = Grids::StructuredCylinderGrid2DAxisymmetric>
        struct Field // : public Grid_t //::CoordinateSystem_t
        {
            static_assert(Grid_t::Dim() == 2ull);
            static_assert(
                std::is_same<
                    typename Grid_t::Axes1,
                    CoordinateTypes::Z>::value);

            Field(
                // propery is along the FIRST Axes1 == Z
                const Logs::StepPropertyGrid &property,
                const Grid_t &grid) noexcept
                : grid{grid},
                  its_values(grid.first_coord.size(), grid.second_coord.size())
            {
                assert(property.log_vals.size() == its_values.rows());
                // extrapolate as const value in the Axes2 direction,
                // though, may be avoided. Element access interface through
                // operator()(i, j) is expected.
                for (std::ptrdiff_t col{0ll}; col < its_values.cols(); ++col)
                    its_values.col(col) = property.log_vals;

                assert(its_values.cols() > 0ll);
                for (std::ptrdiff_t row{0ll}; row < its_values.rows(); ++row)
                    for (std::ptrdiff_t col{1ll}; col < its_values.cols(); ++col)
                        assert(its_values(row, 0) == its_values(row, col));
            }

            const auto& values() const{return its_values;}
            
        protected:
            const Grid_t grid;
            GridNodeValues2D its_values;
        };

        struct Porosity
            : public Field<Grids::StructuredCylinderGrid2DAxisymmetric>
        {
            using Field<Grids::StructuredCylinderGrid2DAxisymmetric>::Field;
        };
        struct Permeability
            : public Field<Grids::StructuredCylinderGrid2DAxisymmetric>
        {
            using Field<Grids::StructuredCylinderGrid2DAxisymmetric>::Field;
        };
        struct HeatVolumetricCapacity
            : public Field<Grids::StructuredCylinderGrid2DAxisymmetric>
        {
            using Field<Grids::StructuredCylinderGrid2DAxisymmetric>::Field;
        };

        template <typename Grid_t = Grids::StructuredCylinderGrid2DAxisymmetric>
        struct FaceInterpolatedField : public Field<Grid_t>
        {
            FaceInterpolatedField(
                const Logs::FaceInterpolatedProperty<typename Grid_t::Axes1> &property,
                const Grid_t &grid)
                : Field<Grid_t>{property, grid},
                  face_vals_axes1(
                      property.face_values.size(),
                      grid.second_coord.size()),
#pragma region AXES2-FACEVALUES
                  face_vals_axes2{
                      Logs::FaceInterpolator<
                          typename Grid_t::Axes2>::interpolate(property,
                                                               grid.second_coord)}
#pragma endregion
            {
#pragma region AXES1-FACEVALUES
                assert(property.face_values.size() >= 0ll);
                // extrapolate Axes1-facevals as const value in the Axes2 direction
                for (std::ptrdiff_t col{0ll}; col < face_vals_axes1.cols(); ++col)
                    face_vals_axes1.col(col) = property.face_values;

                assert(face_vals_axes1.cols() > 0ll);
                for (std::ptrdiff_t row{0ll}; row < face_vals_axes1.rows(); ++row)
                    for (std::ptrdiff_t col{1ll}; col < face_vals_axes1.cols(); ++col)
                        assert(face_vals_axes1(row, 0) == face_vals_axes1(row, col));
#pragma endregion
            }

            // interpolated values at faces normal to Axes1
            FaceValuesContainer face_vals_axes1;
            // interpolated values at faces normal to Axes2
            FaceValuesContainer face_vals_axes2;
        };

        struct HeatConductivity
            : public FaceInterpolatedField<
                  Grids::StructuredCylinderGrid2DAxisymmetric>
        {
            using FaceInterpolatedField<
                Grids::StructuredCylinderGrid2DAxisymmetric>::FaceInterpolatedField;
        };

    } // Properties
} // GPN
