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
    namespace Properties
    {
        using GridNodeValues2D = Eigen::ArrayXX<RealType>;
        using FaceAxes1VauesContainer = GridNodeValues2D;
        using FaceAxes2VauesContainer = GridNodeValues2D;
    } // Properties

    namespace Logs
    {
        template </*CoordinateTypes::R_CylCoord*/>
        struct FaceInterpolator<CoordinateTypes::R_CylCoord>
        {
            using Axes = CoordinateTypes::R_CylCoord;

            const Properties::FaceAxes2VauesContainer face_values;

            template <typename Grid_t>
            FaceInterpolator(
                const StepPropertyGrid &log,
                const Grid_t &grid)
                : face_values{
                      interpolate(log, grid)}
            {
                static_assert(std::is_same<Axes, typename Grid_t::Axes>::value);
            }

        protected:
            template <typename Grid_t>
            static auto interpolate(
                const StepPropertyGrid &log,
                const Grid_t &grid)
            {
                // if dual_size() == 2 --- no internal faces
                assert(grid.dual_size() >= 2ll);

                Properties::FaceAxes2VauesContainer out(log.grid.mesh_size(), grid.dual_size() - 2ll);

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
        template <typename Grid_t = Grids::StructuredCylinderGrid2D>
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
                  values(grid.first_coord.size(), grid.second_coord.size())
            {
                assert(property.log_vals.size() == values.rows());
                // extrapolate as const value in the Axes2 direction,
                // though, may be avoided. Element access interface through
                // operator()(i, j) is expected.
                for (std::ptrdiff_t col{0ll}; col < values.cols(); ++col)
                    values.col(col) = property.log_vals;

                assert(values.cols() > 0ll);
                for (std::ptrdiff_t row{0ll}; row < values.rows(); ++row)
                    for (std::ptrdiff_t col{1ll}; col < values.cols(); ++col)
                        assert(values(row, 0) == values(row, col));
            }

        protected:
            const Grid_t grid;
            GridNodeValues2D values;
        };

        struct Porosity
            : public Field<Grids::StructuredCylinderGrid2D>
        {
            using Field<Grids::StructuredCylinderGrid2D>::Field;
        };
        struct Permeability
            : public Field<Grids::StructuredCylinderGrid2D>
        {
            using Field<Grids::StructuredCylinderGrid2D>::Field;
        };
        struct HeatVolumetricCapacity
            : public Field<Grids::StructuredCylinderGrid2D>
        {
            using Field<Grids::StructuredCylinderGrid2D>::Field;
        };

        template <typename Grid_t = Grids::StructuredCylinderGrid2D>
        struct FaceInterpolatedField : public Field<Grid_t>
        {
            FaceInterpolatedField(
                const Logs::FaceInterpolatedProperty<typename Grid_t::Axes1> &property,
                const Grid_t &grid)
                : Field<Grid_t>{property, grid},
                  axes1_face_vals(
                    property.face_values.size(), 
                    grid.second_coord.size()),
#pragma region AXES2-FACEVALUES
                  axes2_face_vals{
                      Logs::FaceInterpolator<
                          typename Grid_t::Axes2>{
                          property,
                          grid.second_coord}
                          .face_values}
#pragma endregion
            {
#pragma region AXES1-FACEVALUES
                assert(property.face_values.size() >= 0ll);
                // extrapolate Axes1-facevals as const value in the Axes2 direction
                for (std::ptrdiff_t col{0ll}; col < axes1_face_vals.cols(); ++col)
                    axes1_face_vals.col(col) = property.face_values;

                assert(axes1_face_vals.cols() > 0ll);
                for (std::ptrdiff_t row{0ll}; row < axes1_face_vals.rows(); ++row)
                    for (std::ptrdiff_t col{1ll}; col < axes1_face_vals.cols(); ++col)
                        assert(axes1_face_vals(row, 0) == axes1_face_vals(row, col));
#pragma endregion
            }

            FaceAxes1VauesContainer axes1_face_vals;
            FaceAxes2VauesContainer axes2_face_vals;
        };

        struct HeatConductivity
            : public FaceInterpolatedField<
                  Grids::StructuredCylinderGrid2D>
        {
            using FaceInterpolatedField<
                Grids::StructuredCylinderGrid2D>::FaceInterpolatedField;
        };

    } // Properties
} // GPN
