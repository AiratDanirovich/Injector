#pragma once

#include <vector>

#include <Injector/Solver/State1D.hpp>

namespace GPN
{
    namespace Wells
    {

        template <
            typename Grid2D_t,
            typename CrossFlow_t>
        struct SomeWell : public CrossFlow_t
        {
            using Base = CrossFlow_t;

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

            SomeWell(
                const CrossFlow_t &well_base,
                const cptr<Grid2D_t> grid2D_rocks)
                : Base{well_base},
                  grid2D_rocks{grid2D_rocks}
            {
                static_assert(Grid2D_t::l_margin == 3ll);
            }

            template <typename HistoryRecord_t>
            void set_well_flow_field(
                const HistoryRecord_t &record,
                const auto &collector_pressure)
            {
            }

            const cptr<Grid2D_t> grid2D_rocks;
        };
    } // Wells
} // GPN