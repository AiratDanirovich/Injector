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
    namespace Properties
    {
        /// @brief Interpolate physical property between nodes of 2D grid
        /// @tparam Grid_t 2D grid
        template <typename Grid_t>
        struct Field
        {
            static_assert(Grid_t::Dim() == 2ull);
            static_assert(
                std::is_same<
                    typename Grid_t::Axes1,
                    CoordinateTypes::Z>::value);

            Field(const Field &) noexcept = default;
            Field(Field &&) noexcept = default;

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

            auto col(auto i)
            {
                return its_values.col(i);
            }
            const auto &col(auto i) const
            {
                return its_values.col(i);
            }
            
            auto rows() const
            {
                return its_values.rows();
            }

            using Grid_type = Grid_t;

            GridNodeValues2D its_values;
            const cptr<Grid_type> grid;
        };

        template <typename Grid2D_t>
        auto operator/(const Field<Grid2D_t> &lhs, const Field<Grid2D_t> &rhs)
        {
            return Field{lhs.values() / rhs.values(), lhs.grid};
        }

        struct FieldFactory
        {
            template <typename Grid_t>
            static Field<Grid_t> create(
                // property is along the FIRST Axes1 == Z
                const Logs::StepPropertyGrid &property,
                const cptr<Grid_t> grid)
            {
                Field<Grid_t> field{interpolate_r(property, grid), grid};
                // check that interpolation was correct
                for (std::ptrdiff_t col{1ll}; col < field.values().cols(); ++col)
                    for (std::ptrdiff_t row{0ll}; row < field.values().rows(); ++row)
                        assert(field.its_values(row, 0) == field.value(row, col));

                return field;
            }

            template <typename Container_t, typename Grid_t>
            static void apply(
                Field<Grid_t> &field,
                const std::vector<Container_t> &logs)
            {
                for (auto i{0ll}; i < logs.size(); ++i)
                    field.col(i) = logs[i];
                return field;
            }

            template <typename Container_t, typename Grid_t>
            static Field<Grid_t> create(
                const Container_t &base_log, const std::vector<Container_t> &logs,
                const cptr<Grid_t> grid)
            {
                Field temp{create(base_log, grid)};
                apply(temp, logs);
                return temp;
            }

        protected:
            // extrapolate as const value in the Axes2 direction,
            // though, may be avoided. Element access interface through
            // operator()(i, j) is expected.
            template <typename Grid_t>
            static auto interpolate_r(
                const Logs::StepPropertyGrid &property,
                const cptr<Grid_t> grid)
            {
                assert(property.log_vals.size() == grid->first_coord.mesh_size());
                assert(grid->second_coord.mesh_size() > 0ll);

                GridNodeValues2D out(grid->first_coord.mesh_size(), grid->second_coord.mesh_size());
                out.colwise() = property.log_vals;
                return out;
            }
        };

#pragma region ROCKS-PROPERTIES
        template <typename Grid2D_t>
        using Porosity = Field<Grid2D_t>;

        template <typename Grid2D_t>
        using Permeability = Field<Grid2D_t>;
#pragma endregion
#pragma region HEAT-PROPERTIES
        template <typename Grid2D_t>
        using MediumHeatVolumetricCapacity = Field<Grids::StructuredCylinderGrid2DAxisymmetric>;

        template <typename Grid2D_t>
        using HeatConductivity = Field<Grid2D_t>;

        template <typename Grid2D_t>
        struct ThermalDiffusivity
            : public Field<Grid2D_t>
        {
            ThermalDiffusivity(const Properties::MediumHeatVolumetricCapacity<Grid2D_t> &capacity,
                               const Properties::HeatConductivity<Grid2D_t> &conductivity)
                : Field<Grid2D_t>{
                      conductivity / capacity}
            {
            }
        };
#pragma endregion
    } // Properties
} // GPN
