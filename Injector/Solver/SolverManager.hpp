#pragma once

#include <limits>

#include <Injector/History/History.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

namespace GPN
{
    namespace EqSolver
    {
        template <
            typename Grid_t,
            typename Capacity_t>
        struct SolverManager
        {
            using Solver_t = SplittingMethod::Solver<
                Grid_t,
                Capacity_t>;

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
                // loop over the History intries
                for (ptrdiff_t t_step{0ll}; t_step < time_intervals.size(); ++t_step)
                {
                    // calculate the accurate time step
                    size_t internal_step_count{static_cast<size_t>(std::abs(std::ceil(time_intervals[t_step] / numerical_step)))};
                    RealType step{time_intervals[t_step] / internal_step_count};
                    // loop within the histroy entry
                    // with the numerical scheme time step
                    for (size_t id{0ull}; id < internal_step_count; ++id)
                        solver->advance(step);
                }
            }

            const History history;
            const cptr<Solver_t> solver;
        };
    } // EqSolver
} // GPN