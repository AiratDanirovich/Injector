#pragma once

#include <Injector/Grids/Defines.h>

#include <Injector/History/History.hpp>

#include <Injector/Properties/Logs.hpp>

#include <Injector/Model/Phases/PhaseProperties.hpp>

#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
#pragma region INITIAL-CONDITION
        /// @brief Initial pressure is assumed to be constant
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

        template <typename Grid2D_t>
        auto ICFactory(RealType t0, const cptr<Grid2D_t> grid, const auto &ext_pressure)
        {
            return EqSolver::State::State2D{
                EqSolver::State::State2D::FillWithFunctor(
                    *grid, FunctorIC{ext_pressure}, t0)};
        }
#pragma endregion
        template <typename Solver_t>
        struct CompressibleFluidSolver
        {
            CompressibleFluidSolver(
                const ptr<Solver_t> solver)
                : solver{solver}
            {
            }

            void advance(
                const RealType time_step)
            {
                solver->advance(time_step);
            }

            const auto get_state() const
            {
                return solver->get_state();
            }

            const ptr<Solver_t> solver;
        };
    } // Hydrodynamic
} // GPN