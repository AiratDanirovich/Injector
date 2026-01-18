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

        template <typename Grid2D_t>
        auto ICFactory(RealType t0, const cptr<Grid2D_t> grid, const auto &ext_pressure)
        {
            return EqSolver::State::State2D{
                EqSolver::State::State2D::FillWithFunctor(
                    *grid, FunctorIC{ext_pressure}, t0)};
        }
#pragma endregion
#pragma region BOUNDARY-CONDITION
        template <typename Well_t, typename History_t, typename Grid2D_t>
        struct FunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
        {
            //    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
            FunctorBC(
                const cptr<History_t> history,
                const Logs::ExternalPressure &ext_pressure,
                const StepPropertyContainer &rfp,
                const cptr<const Grid2D_t> grid_ptr)
                : history{history},
                  ext_pressure{ext_pressure},
                  grid_ptr{grid_ptr},
                  rfp{rfp}
            {
                assert(rfp.sum() == 1.0);
            }

            RealType operator()(const ptrdiff_t z_id, const RealType r, const RealType,
                                const BCType bc_type) const override
            {
                if (r == grid_ptr->second_coord().dual_front())
                {
                    assert(bc_type == BCType::second);
                    const auto out{rfp(z_id) * history->rate()};
                    return out;
                }

                if (r == grid_ptr->second_coord().dual_back())
                {
                    assert(bc_type == BCType::first);
                    return ext_pressure(z_id);
                }

                assert(false);
                return 0.0;
            }

            RealType operator()(const RealType z, const ptrdiff_t, const RealType,
                                const BCType bc_type) const override
            {
                // boundary conditions are set exactly at domain boundaries
                assert((z == grid_ptr->first_coord().dual_front()) || (z == grid_ptr->first_coord().dual_back()));
                // assume zero diffusion flux in hydrodynamic equation
                assert(bc_type == BCType::second);
                return 0.0;
            }

        protected:
            const Logs::ExternalPressure &ext_pressure;
            const cptr<const Grid2D_t> grid_ptr;
            const StepPropertyContainer rfp;
            const cptr<History_t> history;
        };

        // struct HydroBC : public BoundaryConditions::GeneralBC
        // {
        //     template <typename Grid2D_t>
        //     HydroBC(const cptr<Grid2D_t> &grid,
        //             const cptr<const BCFunctorBase> functor)
        //         : BoundaryConditions::GeneralBC{grid, functor, BoundaryCondition::BCType::second}
        //     {
        //         bc_types[east_id] = BoundaryCondition::BCType::first;
        //     }
        // void set_bc_type(const RealType t)
        // {
        //     this->t = t;
        // }
        // };
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