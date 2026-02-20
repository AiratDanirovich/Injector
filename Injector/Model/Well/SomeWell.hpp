#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>

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
                  grid2D_rocks{grid2D_rocks},
                  flow_axes1_value{
                      FaceValuesContainer::Zero(
                          grid2D_rocks->first_coord().mesh_size() + 1ll,
                          Grid2D_t::l_margin)},
                  flow_axes2_value{
                      FaceValuesContainer::Zero(
                          grid2D_rocks->first_coord().mesh_size(),
                          Grid2D_t::l_margin + 1ll)}
            {
                static_assert(Grid2D_t::l_margin == 3ll);
            }

            void set_well_flow_field(
                const auto &wfp,
                const auto &rfp)
            {

                // set radial flux
                flow_axes2_value.col(0ll) = 0.0;
                flow_axes2_value.col(1ll) = wfp;
                flow_axes2_value.col(2ll) = wfp;
                flow_axes2_value.col(3ll) = rfp;
            }

            FaceValuesContainer flow_axes1_value, flow_axes2_value;
            const cptr<Grid2D_t> grid2D_rocks;
        };
    } // Wells
} // GPN