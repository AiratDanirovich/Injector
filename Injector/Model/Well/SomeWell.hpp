#pragma once

#include <vector>
#include <numeric>
#include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Solver/State1D.hpp>

namespace GPN
{
    namespace Wells
    {
        template <
            typename History_t,
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
                    bottom_pressure.reserve(history->size());
                    bottom_rate.reserve(history->size());

                    rfp.reserve(history->size());
                    wfp.reserve(history->size());
                    verticle_cement_flow.reserve(history->size());
                    verticle_well_flow.reserve(history->size());
                //    cumulative_well_rate.reserve(history->size());
                    well_pressure.reserve(history->size());
                }

                std::vector<RealType> times, bottom_pressure, bottom_rate;
                std::vector<EqSolver::State::State1D>
                    rfp, wfp,
                    verticle_cement_flow,
                    verticle_well_flow,
                //    cumulative_well_rate,
                    well_pressure;
            };

            Solution solution;

            SomeWell(
                const Properties::Rocks::RocksProps<Grid2D_t> &
                    rock_field_props,
                const CrossFlow_t &well_base,
                const cptr<History_t> history,
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
                          Grid2D_t::l_margin + 1ll)},
                  history{history},
                  PI{set_productivity_index(rock_field_props, grid2D_rocks)},
                  rock_field_props{rock_field_props},
                  is_permeable{rock_field_props.base_hydrodynamics.is_permeable},
                  bottom_pressure{0.0},
                  bottom_rate{0.0},
                  solution{history},
                  well_pressure{StepPropertyContainer::Constant(grid2D_rocks->first_coord().mesh_size(), 0.0)}
            {
                static_assert(Grid2D_t::l_margin == 3ll);

                assert(std::all_of(PI.cbegin(), PI.cend(), [](const RealType v)
                                   { return v >= 0.0; }));
                assert(std::any_of(PI.cbegin(), PI.cend(), [](const RealType v)
                                   { return v > 0.0; }));
            }

            void set_well_flow_field(
                double t, RealType t_step,
                const auto &vert_well_flow,
                const auto &vert_cem_flow,
                const auto &wfp,
                const auto &rfp,
                const auto p_bot,
                const auto q_bot,
                const auto well_p)
            {
                const_cast<RealType &>(bottom_pressure) = p_bot;
                const_cast<RealType &>(bottom_rate) = q_bot;
                const_cast<StepPropertyContainer &>(well_pressure) = well_p;

                // set verticle flux
                flow_axes1_value.col(0ll) = vert_well_flow;
                flow_axes1_value.col(1ll) = 0.0;
                flow_axes1_value.col(2ll) = vert_cem_flow;

                // set radial flux
                flow_axes2_value.col(0ll) = 0.0;
                flow_axes2_value.col(1ll) = wfp;
                flow_axes2_value.col(2ll) = wfp;
                flow_axes2_value.col(3ll) = rfp;
#pragma region PUSH-TO-SOLUTION
                {
                    const auto mid_time{history->get_current_record().mid_time};
                    if ((t < mid_time) && (t + t_step > mid_time))
                    {
                        solution.times.push_back(t + t_step / 2.0);
                        solution.bottom_pressure.push_back(bottom_pressure);
                        solution.bottom_rate.push_back(bottom_rate);
                        solution.rfp.emplace_back(RFP());
                        solution.wfp.emplace_back(WFP());
                        solution.verticle_cement_flow.emplace_back(verticle_cement_flow());
                        solution.verticle_well_flow.emplace_back(verticle_well_flow());
                    //    solution.cumulative_well_rate.emplace_back(cumulative_well_rate());
                        solution.well_pressure.emplace_back(well_pressure);
                    }
                }
#pragma endregion
            }

            const RealType z_ref() const
            {
                return history->z_ref;
            }

#pragma region WELL-PARAMS-OF-INTEREST
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

            const StepPropertyContainer cumulative_well_rate() const
            {
                const auto wfp{WFP()};
                StepPropertyContainer out{StepPropertyContainer::Zero(wfp.size() + 1ll)};
                std::partial_sum(
                    wfp.cbegin(), wfp.cend(), out.begin() + 1ll,
                    std::plus<RealType>{});
                assert(out(0ll) == 0.0);
                return bottom_rate - out;
            }
            const RealType bottom_pressure, bottom_rate;

            const StepPropertyContainer well_pressure;
#pragma endregion

            FaceValuesContainer flow_axes1_value, flow_axes2_value;
            const cptr<Grid2D_t> grid2D_rocks;
            const StepPropertyContainer PI;
            const cptr<History_t> history;
            const Properties::Rocks::RocksProps<Grid2D_t> &
                rock_field_props;
            const Logs::IsPermeable &is_permeable;

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

            const auto set_RFP(
                const auto &well_pres_prof,
                const auto &collector_pressure) const
            {
                return Logs::RFPFactory::create_from_container<Logs::RFP>(
                    StepPropertyContainer{(
                                              PI * (well_pres_prof -
                                                    collector_pressure.col(0ll)))
                                              .eval()},
                    is_permeable);
            }

            template <typename HistoryRecord_t>
            const auto &get_RFP(const HistoryRecord_t &record) const
            {
                return Base::rfp;
            }
            template <typename HistoryRecord_t>
            const auto &get_WFP(const HistoryRecord_t &record) const
            {
                return Base::wfp;
            }
            template <typename HistoryRecord_t>
            const auto &get_verticle_cement_flow(const HistoryRecord_t &record) const
            {
                return Base::verticle_flux_in_cement;
            }
            template <typename HistoryRecord_t>
            const auto &get_verticle_well_flow(const HistoryRecord_t &record) const
            {
                return Base::verticle_flux_in_well;
            }
        };
    } // Wells
} // GPN