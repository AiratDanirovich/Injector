#pragma once

#include <numbers>

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

        template <typename Grid2D_t, typename Fluid_t, typename Well_t>
        struct IncompressibleFluidField
        {
            IncompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Logs::ExternalPressure &ext_pressure,
                const Well_t& well,
                const cptr<Grid2D_t> grid2D)
                : fluid{fluid},
                  permeability{permeability},
                  ext_pressure{ext_pressure},
                  P{PressureField::ICFactory(start_time, grid2D, ext_pressure)},
                  thickness_log{grid2D->first_coord.control_volumes}
            {
                const auto& r_grid{grid2D->second_coord};
                const auto r_max{r_grid.dual_back()};
                const auto pi{std::numbers::pi_v<RealType>};

                const StepPropertyContainer temp1{-fluid.viscosity/(2.0*pi)*(r_grid.mesh_nodes/r_max).log()};
                const StepPropertyContainer temp2{thickness_log*permeability.log_vals};
            }

            void set_pressure_field(
                double t, RealType t_step)
            {
            }

        private:
            PressureField P;

            const Fluid_t fluid;
            const Well_t well;
            const cptr<Grid2D_t> grid2D;
            const ControlVolumesContainer &thickness_log;
            const Logs::Permeability permeability;
            const Logs::ExternalPressure ext_pressure;
        };

    } // Hydrodynamic

} // GPN