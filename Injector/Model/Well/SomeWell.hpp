#pragma once

#include <vector>
#include <Eigen/Dense>

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
                const auto &vert_well_flow,
                const auto &vert_cem_flow,
                const auto &wfp,
                const auto &rfp)
            {
                // set verticle flux
                flow_axes1_value.col(0ll) = vert_well_flow;
                flow_axes1_value.col(1ll) = 0.0;
                flow_axes1_value.col(2ll) = vert_cem_flow;

                // set radial flux
                flow_axes2_value.col(0ll) = 0.0;
                flow_axes2_value.col(1ll) = wfp;
                flow_axes2_value.col(2ll) = wfp;
                flow_axes2_value.col(3ll) = rfp;
            }

            const auto RFP() const
            {
                return flow_axes2_value.col(3ll);
            }
            const auto WFP() const
            {
                return flow_axes2_value.col(1ll);
            }
            const auto verticle_cement_flow() const
            {
                return flow_axes1_value.col(2ll);
            }
            const auto verticle_well_flow() const
            {
                return flow_axes1_value.col(0ll);
            }

            FaceValuesContainer flow_axes1_value, flow_axes2_value;
            const cptr<Grid2D_t> grid2D_rocks;

        protected:
            static auto set_productivity_index(
                const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                const ptr<const Grid2D_t> grid2D_rocks)
            {
                // center of cell next to the well
                const auto r3{grid2D_rocks->second_coord().mesh_nodes(0ll)};
                // sandface
                const auto r2_face{grid2D_rocks->second_coord().dual_nodes(0ll)};
                const auto two_pi{2.0 * std::numbers::pi_v<RealType>};
                const auto is_permeable{rock_field_props.base_hydrodynamics.is_permeable.log_vals};

                Eigen::ArrayX<RealType> out{
                    two_pi / std::log(r3 / r2_face) *
                    rock_field_props.mobility_axes2.col(Grid2D_t::l_margin) *
                    grid2D_rocks->first_coord().volumes()};

                for (auto row{0ll}; row < grid2D_rocks->first_coord().mesh_size(); ++row)
                {
                    if (is_permeable(row) == 0.0)
                        out(row) = 0.0;
                }
                return out;
            }
        };
    } // Wells
} // GPN