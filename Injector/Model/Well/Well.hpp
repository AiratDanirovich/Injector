#pragma once

#include <vector>
#include <algorithm>
#include <numbers>
#include <numeric>
#include <cmath>
#include <cassert>

#include <Eigen/Core>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well/WellHoles.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/LogsFactory.hpp>

namespace GPN
{
    struct Friction
    {
    public:
        RealType friction_factor(RealType rate) const noexcept
        {
            return (RealType)0.03;
        }

        RealType f_factor;
    };

    /// @brief Descriptor of materials that fill the
    /// rings that form the well up to sandface.
    /// There is a variation along the verticle direction.
    /// Tube ends at the depth tube_depth.
    struct WellMaterial
    {
        WellMaterial(
            const WellHoles &well_holes,
            const RealType tube_depth,
            const RealType tube_lambda,
            const StationaryPhaseProperties &annulus,
            const StationaryPhaseProperties &cement)
            : well_holes{well_holes},
              tube_lambda{tube_lambda},
              tube_depth{tube_depth},
              annulus{annulus},
              cement{cement}
        {
        }

        /// @brief Depth of the tube
        const RealType tube_depth;
        /// @brief well concentric geometry
        const WellHoles well_holes;
        const RealType tube_lambda;

        const StationaryPhaseProperties annulus, cement;
    };

    struct IWellDesign
    {
        IWellDesign(
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated,
            const StepPropertyContainer &RFP_weights)
            : is_permeable{is_permeable},
              is_perforated{is_perforated},
              RFP_weights{RFP_weights},
              WFP_weights{wfp_weights(Logs::IsGhostLayerFactory::create(is_permeable, is_perforated), RFP_weights, layer_id(is_perforated))},
              weights_sum{RFP_weights.sum()},
              its_top_collector_cell_id{layer_id(is_perforated)}
        {
            assert(is_permeable.size() == is_permeable.grid.dual_steps.size());
            assert(is_perforated.size() == is_perforated.grid.dual_steps.size());
        }
        using Grid_t = Logs::StepPropertyGrid::Grid_t;

        LogValuesContainer get_RFP(const auto &history_record) const = delete;
        LogValuesContainer get_WFP(const auto &history_record) const = delete;
        LogValuesContainer get_verticle_cement_flow(const auto &history_record) const = delete;

        const Logs::IsPermeable is_permeable;
        const Logs::IsPerforated is_perforated;

        const StepPropertyContainer RFP_weights;
        const StepPropertyContainer WFP_weights;
        const RealType weights_sum;

        ptrdiff_t top_collector_cell_id() const
        {
            return its_top_collector_cell_id;
        }

    private:
        static StepPropertyContainer wfp_weights(
            const Logs::IsGhostLayer &is_ghost,
            const StepPropertyContainer &RFP_weights,
            const ptrdiff_t top_collector_cell_id)
        {
            // the well rate is zero at the ghost layer
            StepPropertyContainer WFP_weights{RFP_weights * (1.0 - is_ghost.log_vals)};
            // all ghost layer fluxes flow through the top collector layer
            WFP_weights(top_collector_cell_id) = (RFP_weights * is_ghost.log_vals).sum() + RFP_weights(top_collector_cell_id);
            return WFP_weights;
        }

        static ptrdiff_t layer_id(const auto &indicator)
        {
            const auto perforated_it = std::ranges::find(indicator.log_vals, 1.0);
            return std::distance(indicator.log_vals.cbegin(), perforated_it);
        }
        // cell id where the cross-flow leaves the column
        const ptrdiff_t its_top_collector_cell_id{-1ll};
    };

    struct Well_Explicit
        : public IWellDesign
    {
        Well_Explicit(
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated,
            const StepPropertyContainer &weights)
            : IWellDesign{
                  is_permeable,
                  is_perforated,
                  /*RFP_weights*/ weights}
        {
        }

        template <typename HistoryRecord_t>
        auto get_RFP(const HistoryRecord_t &history_record) const
        {
            return get_RFP(history_record.rate, history_record.pressure);
        }
        template <typename HistoryRecord_t>
        auto get_WFP(const HistoryRecord_t &history_record) const
        {
            return get_WFP(history_record.rate, history_record.pressure);
        }
        template <typename HistoryRecord_t>
        auto get_verticle_cement_flow(const HistoryRecord_t &history_record) const
        {
            return get_verticle_cement_flow(history_record.rate, history_record.pressure);
        }

    protected:
        StepPropertyContainer get_verticle_cement_flow(
            RealType rate,
            RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                throw std::invalid_argument("cement flow: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                //        assert(rate >= 0.0);

                const auto rfp_vals{get_RFP(rate, pressure)};
                // leftover flowrate along the well
                LogValuesContainer cum_sum{LogValuesContainer::Zero(is_permeable.grid.dual_size())};
                std::partial_sum(
                    rfp_vals.cbegin(),
                    rfp_vals.cbegin() + top_collector_cell_id(),
                    cum_sum.begin() + 1ll, std::plus<RealType>());
                // the flow in cement goes UP,
                // while the flow in the tube -- DOWN
                return -cum_sum;
            }
            else
                throw std::invalid_argument("cement flow: Either rate or pressure must be set, but not both.");
        }

        StepPropertyContainer get_RFP(
            RealType rate,
            RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                throw std::invalid_argument("RFP: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                //    assert(rate >= 0.0);
                return ((rate / weights_sum) * RFP_weights).eval();
            }
            else
                throw std::invalid_argument("RFP: Either rate or pressure must be set, but not both.");
        }

        StepPropertyContainer get_WFP(
            RealType rate,
            RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                throw std::invalid_argument("WFP: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                //    assert(rate >= 0.0);
                return ((rate / weights_sum) * WFP_weights).eval();
            }
            else
                throw std::invalid_argument("WFP: Either rate or pressure must be set, but not both.");
        }
    };

    /// @brief Assumption: multiple flows can start from the same well node,
    /// but they have to come to different layers of the reservoir.
    /// This allows to assemble the WFP from the RFP
    /// @param is_permeable
    /// @param is_perforated
    /// @param RFP_weights
    /// @param cross_flows
    /// @return
    Logs::WFP_weights create_WFP_weights(
        const Logs::IsPerforated &is_perforated,
        const Logs::RFP_weights &RFP_w,
        const CrossFlow::CrossFlows &cross_flows)
    {
        StepPropertyContainer wfp_step_prop_grid{(is_perforated * RFP_w).log_vals};
        const auto is_damaged{Logs::IsDamagedFactory::create(
            cross_flows,
            is_perforated.log_vals,
            is_perforated.grid)};

        for (const auto &cf : cross_flows.cross_flow_data)
            wfp_step_prop_grid(cf.from_id) += RFP_w(cf.to_id);

        assert(wfp_step_prop_grid.sum() == RFP_w.log_vals.sum());

        assert(is_perforated.size() == wfp_step_prop_grid.size());
        for (auto i{0ll}; i < wfp_step_prop_grid.size(); ++i)
        {
            assert(
                ((is_perforated(i) != is_damaged(i)) &&
                 (wfp_step_prop_grid(i) > 0.0)) ||
                ((is_perforated(i) == 0.0) && (is_damaged(i) == 0.0) &&
                 (wfp_step_prop_grid(i) == 0.0)));
        }

        return Logs::WFPFactory::create_from_container(
            wfp_step_prop_grid, is_perforated + is_damaged);
    }

    struct Well_CrossFlow
    {
        Well_CrossFlow(
            const Logs::RFP &RFP_w,
            const Logs::WFP &WFP_w,
            const CrossFlow::CrossFlows &cross_flows)
            : RFP_weights{RFP_w.log_vals},
              weights_sum{RFP_w.log_vals.sum()},
              WFP_weights{WFP_w.log_vals},
              cross_flows{cross_flows}
        {
            assert(RFP_weights.log_vals.sum() == WFP_weights.log_vals.sum());
        }

        Well_CrossFlow(
            const Logs::RFP &RFP_w,
            const Logs::WFP &WFP_w,
            const std::vector<RealType> &from_coords,
            const std::vector<ptrdiff_t> &to_layers)
            : Well_CrossFlow{
                  RFP_w, WFP_w,
                  CrossFlow::CrossFlows{RFP_w, from_coords, to_layers}}
        {
        }

        template <typename HistoryRecord_t>
        auto get_RFP(const HistoryRecord_t &history_record) const
        {
            return get_RFP(history_record.rate, history_record.pressure);
        }
        template <typename HistoryRecord_t>
        auto get_WFP(const HistoryRecord_t &history_record) const
        {
            return get_WFP(history_record.rate, history_record.pressure);
        }
        template <typename HistoryRecord_t>
        auto get_verticle_cement_flow(const HistoryRecord_t &history_record) const
        {
            return get_verticle_cement_flow(history_record.rate, history_record.pressure);
        }
        template <typename HistoryRecord_t>
        auto get_verticle_well_flow(const HistoryRecord_t &history_record) const
        {
            return get_verticle_well_flow(history_record.rate, history_record.pressure);
        }

        const StepPropertyContainer &RFP_weights;
        const StepPropertyContainer &WFP_weights;

    protected:
        StepPropertyContainer get_RFP(
            const RealType rate,
            const RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                throw std::invalid_argument("RFP: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                //                assert(rate >= 0.0);
                return ((rate / weights_sum) * RFP_weights).eval();
            }
            else
                throw std::invalid_argument("RFP: Either rate or pressure must be set, but not both.");
        }

        StepPropertyContainer get_WFP(
            const RealType rate,
            const RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                throw std::invalid_argument("WFP: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                //    assert(rate >= 0.0);
                return ((rate / weights_sum) * WFP_weights).eval();
            }
            else
                throw std::invalid_argument("WFP: Either rate or pressure must be set, but not both.");
        }

        StepPropertyContainer get_verticle_cement_flow(
            const RealType rate,
            const RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                throw std::invalid_argument("RFP: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                //                assert(rate >= 0.0);
                return ((rate / weights_sum) * cross_flows.get_normalized_verticle_flux()).eval();
            }
            else
                throw std::invalid_argument("RFP: Either rate or pressure must be set, but not both.");
        }

        StepPropertyContainer get_verticle_well_flow(
            const RealType rate,
            const RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                throw std::invalid_argument("RFP: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));

                const auto wfp{get_WFP(rate, pressure)};
                StepPropertyContainer out(StepPropertyContainer::Zero(wfp.rows() + 1ll));
                std::partial_sum(wfp.cbegin(), wfp.cend(), out.begin() + 1ll, std::plus<RealType>{});
                out = rate - out;

                return out;
            }
            else
                throw std::invalid_argument("verticle well flow: Either rate or pressure must be set, but not both.");
        }

        const RealType weights_sum;

        const CrossFlow::CrossFlows cross_flows;
    };

    /**
     * @brief Incompressible fluid is assumed in reservoir
     */
    struct Well_KH
        : public Well_Explicit
    {
        Well_KH(
            const PhaseProperties &fluid,
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated,
            const StepPropertyContainer &permeability,
            const WellHoles &holes,
            const RealType Rext)
            : Well_Explicit{
                  is_permeable,
                  is_perforated,
                  /*RFP_weights*/ permeability * is_permeable.grid.dual_steps * (StepPropertyContainer)is_permeable},
              fluid{fluid}, log_dist{std::log(Rext / holes.sandface_radius)}
        {
            assert(permeability.size() == is_permeable.grid.dual_steps.size());
        }

        template <typename HistoryRecord_t>
        auto get_RFP(const HistoryRecord_t &history_record) const
        {
            return get_RFP(history_record.rate, history_record.pressure);
        }
        template <typename HistoryRecord_t>
        auto get_WFP(const HistoryRecord_t &history_record) const
        {
            return get_WFP(history_record.rate, history_record.pressure);
        }

    protected:
        StepPropertyContainer get_RFP(
            RealType rate,
            RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                assert(!std::isnan(pressure));
                rate = 2 * std::numbers::pi / fluid.viscosity / log_dist * pressure * weights_sum;
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                pressure = rate / (2 * std::numbers::pi / fluid.viscosity / log_dist * weights_sum);
            }
            else
                assert("Incorrect injector regime!");

            assert(((pressure == 0.0) && (rate == 0.0)) || ((pressure > 0.0) && (rate > 0.0)));
            return ((rate / weights_sum) * RFP_weights).eval();
        }

        StepPropertyContainer get_WFP(
            RealType rate,
            RealType pressure) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                assert(!std::isnan(pressure));
                rate = 2 * std::numbers::pi / fluid.viscosity / log_dist * pressure * weights_sum;
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                pressure = rate / (2 * std::numbers::pi / fluid.viscosity / log_dist * weights_sum);
            }
            else
                assert("Incorrect injector regime!");

            assert(((pressure == 0.0) && (rate == 0.0)) || ((pressure > 0.0) && (rate > 0.0)));

            return ((rate / weights_sum) * WFP_weights).eval();
        }

    private:
        const RealType log_dist;
        const PhaseProperties fluid;
    };
} // GPN
