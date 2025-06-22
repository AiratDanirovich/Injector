#pragma once

#include <vector>
#include <algorithm>
#include <numbers>
#include <cmath>
#include <cassert>

#include <Eigen/Core>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>
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

    struct WellHoles
    {
        WellHoles(
            const RealType tube_radius,
            const RealType sandface_radius)
            : tube_radius{tube_radius},
              sandface_radius{sandface_radius}
        {
            assert(tube_radius < sandface_radius);
        }

        std::vector<RealType> generate_uniform_radial_grid(const RealType r_min, const RealType r_max, const ptrdiff_t r_nodes) const
        {
            assert(r_min < tube_radius);
            assert(r_max > sandface_radius);
            assert(r_nodes > 1ll);

            std::vector<RealType> out;
            out.reserve(r_nodes + 2ll);

            out.push_back(r_min);
            out.push_back(tube_radius);
            out.push_back(sandface_radius);

            const RealType step{(r_max - sandface_radius) / (r_nodes - 1ll)};
            for (auto i{2ll}; i < r_nodes; ++i)
                out.push_back(out.back() + step);
            out.push_back(r_max);

            assert(out.front() == r_min);
            for (auto i{1ull}; i < out.size(); ++i)
                assert(out[i] > out[i - 1ull]);
            assert(out.back() == r_max);

            return out;
        }

        std::vector<RealType> generate_log_radial_grid(
            const RealType r_min, const RealType r_max,
            const RealType q,        // ratio of adjascent steps
            const RealType max_step) // max allowed step
            const
        {
            assert(r_min < tube_radius);
            assert(r_max > sandface_radius);
            assert(max_step > tube_radius);
            assert(max_step > sandface_radius - tube_radius);
            assert(q >= 1.0);
            // base, minimum step for geometric progression
            const RealType base_step = sandface_radius - tube_radius;

            if (q == 1.0)
            {
                // uniform grid
                return generate_uniform_radial_grid(r_min, r_max, (ptrdiff_t)std::ceil((r_max - r_min) / base_step));
            }
            else
            {
                // non-uniform grid
                // nmbr of steps within segment [sandface_radius; r_max],
                // which increase geometrically
                ptrdiff_t nx{
                    (ptrdiff_t)std::ceil(
                        std::log(1.0 + (r_max - sandface_radius) / base_step * (q - 1.0)) /
                        std::log(q))};

                std::vector<RealType> out;
                out.reserve(nx + 20ll);

                out.push_back(r_min);           // push leftmost boundary
                out.push_back(tube_radius);     // push tube radius
                out.push_back(sandface_radius); // push sandface radius

                // recalculate the base step
                const RealType hx{base_step}; //{(r_max - sandface_radius) * (q - 1.0) / (std::pow(q, nx) - 1.0)};
                assert(hx <= base_step);

                if (max_step < hx * (std::pow(q, nx - 1ll)))
                {
                    // the furthest steps are too large
                    for (auto i{0ll}; i < nx; ++i)
                        out.push_back(out.back() + std::min(max_step, hx * std::pow(q, i)));
                    // tail of the segment where the geometric steps are too large
                    while (out.back() < r_max)
                        out.push_back(out.back() + max_step);
                }
                else
                {
                    // the furthest steps are fine
                    for (auto i{0ll}; i < nx; ++i)
                        out.push_back(out.back() + hx * std::pow(q, i));
                }
                for (auto i{1ull}; i < out.size(); ++i)
                    assert(out[i] > out[i - 1ull]);

                return out;
            }
        }

        const RealType tube_radius;
        const RealType sandface_radius;
    };

    struct IWellDesign
    {
        IWellDesign(
            const PhaseProperties &fluid,
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated,
            const StepPropertyContainer &RFP_weights)
            : is_permeable{is_permeable},
              is_perforated{is_perforated},
              RFP_weights{RFP_weights},
              WFP_weights{wfp_weights(Logs::IsGhostLayerFactory::create(is_permeable, is_perforated), RFP_weights, layer_id(is_perforated))},
              fluid{fluid},
              weights_sum{RFP_weights.sum()},
              its_top_collector_cell_id{layer_id(is_perforated)}
        {
            assert(is_permeable.size() == is_permeable.grid.dual_steps.size());
            assert(is_perforated.size() == is_perforated.grid.dual_steps.size());
        }
        using Grid_t = Logs::StepPropertyGrid::Grid_t;

        LogValuesContainer get_RFP(const auto &history_record) const = delete;
        LogValuesContainer get_WFP(const auto &history_record) const = delete;

        const Logs::IsPermeable is_permeable;
        const Logs::IsPerforated is_perforated;
        //    const Logs::IsGhostLayer is_ghost;

        const StepPropertyContainer RFP_weights;
        const StepPropertyContainer WFP_weights;
        const RealType weights_sum;
        const PhaseProperties fluid;

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

    private:
        static ptrdiff_t layer_id(const auto &indicator)
        {
            const auto perforated_it = std::ranges::find(indicator.log_vals, 1.0);
            return std::distance(indicator.log_vals.cbegin(), perforated_it);
        }
        const ptrdiff_t its_top_collector_cell_id{-1ll};
    };

    struct Well_KH
        : public IWellDesign
    {
        Well_KH(
            //    const RealType tube_depth,
            const PhaseProperties &fluid,
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated,
            const StepPropertyContainer &permeability,
            const WellHoles &holes,
            const RealType Rext)
            : IWellDesign{
                  fluid,
                  is_permeable,
                  is_perforated,
                  /*RFP_weights*/ permeability * is_permeable.grid.dual_steps * (StepPropertyContainer)is_permeable},
              log_dist{std::log(Rext / holes.sandface_radius)}
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
    };

    struct Well_Explicit
        : public IWellDesign
    {
        Well_Explicit(
            //    const RealType tube_depth,
            const PhaseProperties &fluid,
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated,
            const StepPropertyContainer &permeability,
            const WellHoles &holes,
            const RealType Rext)
            : IWellDesign{
                  fluid,
                  is_permeable,
                  is_perforated,
                  /*RFP_weights*/ permeability * is_permeable.grid.dual_steps * (StepPropertyContainer)is_permeable}
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
                throw std::invalid_argument("RFP: Rate must be set");
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                assert(rate >= 0.0);
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
                assert(rate >= 0.0);
                return ((rate / weights_sum) * WFP_weights).eval();
            }
            else
                throw std::invalid_argument("WFP: Either rate or pressure must be set, but not both.");
        }
    };
} // GPN
