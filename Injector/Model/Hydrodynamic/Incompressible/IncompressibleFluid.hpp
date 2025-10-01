#pragma once

#include <Injector/Solver/State2D.hpp>

#include <Injector/Model/Phases/PhaseProperties.hpp>



namespace GPN
{
    namespace Hydrodynamic
    {
        struct PressureField 
        : public EqSolver::State::State2D
        {

        };

        struct IncompressibleFluid
        {
            IncompressibleFluid(
                const Water& fluid)
            : fluid{fluid} //, history{}
            {

            }


            private:
            Water fluid;
        };


    } // Hydrodynamic

} // GPN