#pragma once

#include <Injector/Solver/State2D.hpp>

namespace EqSolver
{
    namespace Problem
    {
        struct InitialCondition : public State::State2D
        {
            InitialCondition(
                const State::State2D &state) noexcept
                : State::State2D{state}
            {
            }
        };
    } // Problem
} // EqSolver