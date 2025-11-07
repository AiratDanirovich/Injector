#pragma once

#include <limits>

namespace GPN
{
    namespace EqSolver
    {
        template <
            typename History_t,
            typename Solver_t>
        struct SolverManager
        {
            SolverManager(
                const ptr<History_t> history,
                const ptr<Solver_t> solver)
                : history{history},
                  solver{solver}
            {
            }

            void run(const RealType numerical_step = std::numeric_limits<RealType>::max())
            {
                const ptrdiff_t size{history->size()};
                const auto &time_intervals{history->time_steps};
                // loop over the History entries
                for (ptrdiff_t t_step{0ll}; t_step < size; ++t_step)
                {
                    // set history pointer at the current record,
                    // so it pointed at the last record when the simulation is finished
                    history->advance();
                    // calculate the accurate time step
                    size_t internal_step_count{
                        static_cast<size_t>(
                            std::abs(std::ceil(time_intervals[t_step] / numerical_step)))};
                    if(internal_step_count % 2 == 0)
                        ++internal_step_count;
                    const RealType step{time_intervals[t_step] / internal_step_count};
                    // loop within the histroy entry
                    // with the numerical scheme time step
                    for (size_t id{0ull}; id < internal_step_count; ++id)
                        solver->advance(step);
                    solver->save_state();
                }
                assert(history->record_id() == history->size()-1ll);
            }

            const ptr<History_t> history;
            const ptr<Solver_t> solver;
        };
    } // EqSolver
} // GPN