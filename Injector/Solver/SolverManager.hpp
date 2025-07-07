#pragma once

#include <limits>

#include <Injector/History/History.hpp>

namespace GPN
{
    namespace EqSolver
    {
        template <
            typename Solver_t>
        struct SolverManager
        {
            SolverManager(
                const History &history,
                cptr<Solver_t> solver)
                : history{history},
                  solver{solver}
            {
            }

            void run(RealType numerical_step = std::numeric_limits<RealType>::max())
            {
                const auto &time_intervals{history.time_steps};
                // loop over the History entries
                for (ptrdiff_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
                {
                    // calculate the accurate time step
                    size_t internal_step_count{
                        static_cast<size_t>(
                            std::abs(std::ceil(time_intervals[t_step] / numerical_step)))};
                    RealType step{time_intervals[t_step] / internal_step_count};
                    // loop within the histroy entry
                    // with the numerical scheme time step
                    for (size_t id{0ull}; id < internal_step_count; ++id)
                    {
                        solver->advance(step);
                    }
                    solver->save_state();
                }
            }

            const History history;
            const cptr<Solver_t> solver;
        };
    } // EqSolver
} // GPN