#pragma once

#include <Injector/Solver/State2D.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace Problem
        {
            struct InitialCondition : public State::State2D
            {
                using State::State2D::State2D;
                InitialCondition(
                    const State::State2D &state) noexcept
                    : State::State2D{state}
                {
                }
            };
        } // Problem
    } // EqSolver
} // GPN