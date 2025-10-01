#pragma once

#include <Injector/Solver/State2D.hpp>

#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
        struct PressureField
            : public EqSolver::State::State2D
        {
            PressureField(const EqSolver::State::State2D &state)
                : EqSolver::State::State2D{state}
            {
            }

            /// @brief Initial pressure is assumed to be constant with r
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
            static auto ICFactory(RealType t0, const Grid_t_ptr grid, const Logs::ExternalPressure &ext_pressure)
            {
                return EqSolver::State::State2D{EqSolver::State::State2D::FillWithFunctor(*grid, FunctorIC{ext_pressure}, t0)};
            }
        };

        template <typename Grid2D_t, typename Fluid_t>
        struct IncompressibleFluidField
        {
            IncompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Logs::ExternalPressure &ext_pressure,
                const cptr<Grid2D_t> grid2D)
                : fluid{fluid},
                  permeability{permeability},
                  ext_pressure{ext_pressure},
                  P{PressureField::ICFactory(start_time, grid2D, ext_pressure)}
            {
            }

        private:
            PressureField P;

            const Fluid_t fluid;
            const cptr<Grid2D_t> grid2D;
            const Logs::Permeability permeability;
            const Logs::ExternalPressure ext_pressure;
        };

    } // Hydrodynamic

} // GPN