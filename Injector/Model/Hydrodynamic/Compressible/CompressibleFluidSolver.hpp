#pragma once

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/Logs.hpp>

#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
#pragma region INITIAL-CONDITION
        /// @brief Initial temperature is assumed to be constant
        struct FunctorIC : public GPN::InitialConditions::ICFunctorBase
        {
            FunctorIC(const Logs::ExternalPressure &ext_pressure)
                : ext_pressure{ext_pressure}
            {
            }
            RealType operator()(const ptrdiff_t z_id, const ptrdiff_t, GPN::RealType) const override
            {
                return ext_pressure(z_id);
            }

        protected:
            const Logs::ExternalPressure &ext_pressure;
        };

        template <typename Grid_t_ptr>
        auto ICFactory(RealType t0, const Grid_t_ptr grid, const Logs::Geotherma &geotherma)
        {
            return EqSolver::State::State2D{
                EqSolver::State::State2D::FillWithFunctor(
                    *grid, FunctorIC{geotherma}, t0)};
        }
#pragma endregion
#pragma region BOUNDARY-CONDITION
        template <typename Well_t, typename Hydro_t>
        struct FunctorBC : public BoundaryConditions::BCFunctorBase
        {
            using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
            using ConvectionFieldFactory_t =
                FaceProperties::IncompressibleRatesFactory<
                    Grid2D_t, Well_t, PhasePropertiesJT, Hydro_t>;
            FunctorBC(
                const Logs::ExternalPressure &ext_pressure const cptr<const Grid2D_t> grid_ptr)
                : flow_field{flow_field},
                  ext_pressure{ext_pressure},
                  grid_ptr{grid_ptr}
            {
            }

            RealType operator()(const ptrdiff_t z_id, const RealType r, const RealType t) const override
            {
                if (r == grid_ptr->second_coord.dual_front())
                    return flow_field.get_heat_flow_in_axes2()(z_id, 0ll) * flow_field.get_temperature();

                if (r == grid_ptr->second_coord.dual_back())
                    return 0.0;

                assert(false);
                return 0.0;
            }

            RealType operator()(const RealType z, const ptrdiff_t r_id, const RealType t) const override
            {
                if (z == grid_ptr->first_coord.dual_front())
                    return flow_field.get_heat_flow_in_axes1()(0ll, r_id) * flow_field.get_temperature();

                if (z == grid_ptr->first_coord.dual_back())
                    return 0.0;

                assert(false);
                return 0.0;
            }

        protected:
            const Logs::ExternalPressure ext_pressure;
            const cptr<const Grid2D_t> grid_ptr;
            const ConvectionFieldFactory_t &flow_field;
        };
#pragma endregion
        struct CompressibleFluidSolver
        {
            CompressibleFluidSolver(
                const Logs::ExternalPressure &ext_pressure)
                : initial_condition{FunctorIC{ext_pressure}}
            {
            }

            const FunctorIC initial_condition;
        };
    } // Hydrodynamic

} // GPN