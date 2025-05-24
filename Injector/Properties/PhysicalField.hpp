#pragma once

#include <memory>
#include <type_traits>
#include <omp.h>
#include <cassert>

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
                        static_cast<StepPropertyContainer>(log.log_vals));
                return out;
            }
        };
    } // Logs

    namespace Properties
    {
        /// @brief Interpolate physical property between nodes of 2D grid
        /// @tparam Grid_t 2D grid
        template <typename Grid_t = Grids::StructuredCylinderGrid2DAxisymmetric>
        struct Field
        {
            static_assert(Grid_t::Dim() == 2ull);
            static_assert(
                std::is_same<
                    typename Grid_t::Axes1,
                    CoordinateTypes::Z>::value);

            Field(const Field &) noexcept = default;
            Field(
                // property is along the FIRST Axes1 == Z
                const Logs::StepPropertyGrid &property,
                const cptr<Grid_t> grid) noexcept
                : Field{interpolate_r(property, grid), grid}
            {
                // check that interpolation was correct
                for (std::ptrdiff_t row{0ll}; row < its_values.rows(); ++row)
                    for (std::ptrdiff_t col{1ll}; col < its_values.cols(); ++col)
                        assert(its_values(row, 0) == its_values(row, col));
            }

            Field(const GridNodeValues2D &vals,
                  const cptr<Grid_t> grid)
                : its_values{vals},
                  grid{grid}
            {
                assert(its_values.rows() == grid->first_coord.mesh_size());
                assert(its_values.cols() == grid->second_coord.mesh_size());
                assert(grid->second_coord.mesh_size() > 0ll);
            }

            const auto &values() const { return its_values; }
            RealType value(auto first, auto second) const { return values()(first, second); }

            auto operator/(const Field &rhs) const
            {
                const auto &lhs{*this};
                return Field{lhs.values() / rhs.values(), lhs.grid};
            }

            using Grid_type = Grid_t;

            const GridNodeValues2D its_values;
            const cptr<Grid_type> grid;

            template <typename Container_t>
            static Field set_from_multiple(
                const Container_t &base_log, const std::vector<Container_t> &logs,
                const cptr<Grid_t> &grid)
            {
                auto temp{interpolate_r(base_log, grid)};
                for (auto i{0ll}; i < logs.size(); ++i)
                    temp.col(i) = logs[i];
                return temp;
            }

        private:
            // extrapolate as const value in the Axes2 direction,
            // though, may be avoided. Element access interface through
            // operator()(i, j) is expected.
            static auto interpolate_r(const Logs::StepPropertyGrid &property, const cptr<Grid_t> grid)
            {
                assert(property.log_vals.size() == grid->first_coord.mesh_size());
                assert(grid->second_coord.mesh_size() > 0ll);

                GridNodeValues2D out(grid->first_coord.mesh_size(), grid->second_coord.mesh_size());
                out.colwise() = property.log_vals;
                return out;
            }
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
            // Permeability(
            //     const Logs::StepPropertyGrid &property,
            //     const cptr<Grids::StructuredCylinderGrid2DAxisymmetric> grid)
            //     : Field{property, grid}
            // {
            // }
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
                const cptr<Grid_t> grid)
                : Field<Grid_t>{property, grid},
                  face_vals_axes1(
                      property.face_values.size(),
                      grid->second_coord.mesh_size()),
                  // #pragma region AXES2-FACEVALUES
                  face_vals_axes2{
                      Logs::FaceInterpolator<
                          typename Grid_t::Axes2>::interpolate(property,
                                                               grid->second_coord)}
            // #pragma endregion
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

            using typename Field<Grid_t>::Grid_type;

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

        struct ThermalDiffusivity
            : public Field<
                  Grids::StructuredCylinderGrid2DAxisymmetric>
        {
            ThermalDiffusivity(const Properties::HeatVolumetricCapacity &capacity,
                               const Properties::HeatConductivity &conductivity)
                : Field{
                      conductivity / capacity}
            {
            }
        };

    } // Properties
} // GPN
