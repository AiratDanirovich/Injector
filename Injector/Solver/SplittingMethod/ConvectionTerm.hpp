#pragma once

#include <memory>

#include <Injector/Grids/Grids.hpp>
#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Model/HydrodynamicSolver.hpp>

namespace GPN
{
    namespace EqSolver
    {
    //    /// @brief
    //    /// @tparam PtrHydrodynamic_t Read-only pointer (shared) to evaluator of hydrodynamic
    //    template <typename Grid_t, typename Hydrodynamic_t>
        struct ConvectionTerm
        {
            using Grid_t = Grids::StructuredCylinderGrid2DAxisymmetric;
            using Hydrodynamic_t = Model::Injector::Hydrodynamic2D<Water>;


            ConvectionTerm(
                const Hydrodynamic_t& hydro,
                const Grid_t& grid) noexcept
                : hydrodynamic{hydro}
                , grid{grid}
            {
            }

        protected:
            const Hydrodynamic_t& hydrodynamic;
            const Grid_t& grid;
            // FluxField flux_field;
            // HeatCapacity
        };


















        /// @brief Manages convection term factors at fixed stripe along y-axis. It stores c*rho*volume factor,
        /// and multiplies it by the flow field at cell boundaries on request.
        /// @tparam PtrHydrodynamic_t Read-only pointer (shared) to evaluator of hydrodynamic
        template <typename PtrHydrodynamic_t, typename Grid_t>
        struct ConvectionTermX : public ConvectionTerm<PtrHydrodynamic_t>
        {
            ConvectionTerm(
                PtrHydrodynamic_t hydro,
                const Grid_t &grid) noexcept
                : ConvectionTermX{hydro, grid}
            {
            }

        protected:
            const PtrHydrodynamic_t hydrodynamic;
            FluxField flux_field;

            Conductivity_f conductivity_y_bounds{
                set_volumetricheatcapasity_at_y_bounds(properties->conductivity_f, grid)};

            template <typename Grid_t>
            auto set_volumetricheatcapasity_at_y_bounds(
                const Conductivity_f &conductivity_nodes,
                const Grid_t &grid) const
            {
                Conductivity_f conductivity_y_bounds{
                    grid.first_coord.size(),
                    grid.second_coord.size() - 1};

#pragma omp parallel for
                for (ptrdiff_t i = 0; i < conductivity_y_bounds.rows(); ++i)
                    for (ptrdiff_t j = 0; j < conductivity_y_bounds.cols(); ++j)
                        conductivity_y_bounds(i, j) =
                            (conductivity_nodes(i, j + 1) + conductivity_nodes(i, j)) / 2.0;
                return conductivity_y_bounds;
            }
        };
    } // EqSolver
} // GPN
