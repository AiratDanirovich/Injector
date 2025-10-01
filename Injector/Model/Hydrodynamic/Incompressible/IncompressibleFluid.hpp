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
        };

        template <typename Grid2D_t, typename Fluid_t>
        struct IncompressibleFluidField
        {
            IncompressibleFluidField(
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Logs::ExternalPressure &ext_pressure,
                const cptr<Grid2D_t> grid2D)
                : fluid{fluid},
                  permeability{permeability},
                  ext_pressure{ext_pressure}
            {
            }

        private:
            //            PressureField P;

            const Fluid_t fluid;
            const cptr<Grid2D_t> grid2D;
            const Logs::Permeability permeability;
            const Logs::ExternalPressure ext_pressure;
        };

    } // Hydrodynamic

} // GPN