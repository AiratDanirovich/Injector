#pragma once

#include <vector>

#include <Injector/Solver/State1D.hpp>

namespace GPN
{
    namespace Wells
    {

        struct SomeWell
        {
            struct Solution
            {
                Solution(const auto history)
                {
                    times.reserve(history->size());
                    states.reserve(history->size());
                }

                std::vector<RealType> times;
                std::vector<EqSolver::State::State1D> states;
            };


        };
    } // Wells
} // GPN