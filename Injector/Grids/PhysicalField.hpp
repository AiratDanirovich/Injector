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
                for (std::ptrdiff_t col{0ll}; col < values.cols(); ++col)
                    values.col(col) = property.log_vals;
            }

        protected:
            const Grid_t &grid;
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

        template <typename Grid_t = Grids::StructuredCylinderGrid2D>
        struct FaceInterpolatedField : public Field<Grid_t>
        {
            FaceInterpolatedField(
                const Logs::FaceInterpolatedProperty<typename Grid_t::Axes1> &property,
                const Grid_t &grid)
                : Field<Grid_t>{property, grid},
                  axes1_vals(property.face_values.size(), grid.second_coord.size()),
                  axes2_vals(property.log_vals.size(), grid.second_coord.dual_size() - 2ll)
            {
            }

            FaceAxes1VauesContainer axes1_vals;
            FaceAxes2VauesContainer axes2_vals;
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
