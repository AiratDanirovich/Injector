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

        std::vector<RealType> generate_uniform_radial_grid(const RealType r_min, const RealType r_max, const ptrdiff_t r_nodes)
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
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated)
            : is_permeable{is_permeable},
              is_perforated{is_perforated}
        {
        }
        using Grid_t = Logs::StepPropertyGrid::Grid_t;
        virtual LogValuesContainer get_RFP(RealType rate, RealType pressure) const = 0;

        const Logs::IsPermeable is_permeable;
        const Logs::IsPerforated is_perforated;
    };

    struct Well_KH_FixedRate
        : public IWellDesign
    {
        Well_KH_FixedRate(
            //    const RealType tube_depth,
            const PhaseProperties &fluid,
            const Logs::IsPermeable &is_permeable,
            const Logs::IsPerforated &is_perforated,
            const StepPropertyContainer &permeability,
            const WellHoles &holes,
            const RealType Rext)
            : IWellDesign{is_permeable, is_perforated},
              RFP_weights{permeability * is_permeable.grid.dual_steps * (StepPropertyContainer)is_permeable},
              top_collector_cell_id{layer_id(is_perforated)},
              ghost_layer_cell_id{layer_id(is_permeable)},
              fluid{fluid},
              log_dist{std::log(Rext / holes.sandface_radius)}
        {
            assert(permeability.size() == is_permeable.grid.dual_steps.size());
            assert(is_permeable.size() == is_permeable.grid.dual_steps.size());
            assert(is_perforated.size() == is_perforated.grid.dual_steps.size());

            weights_sum = RFP_weights.sum();
            WFP_weights = RFP_weights;
            // the well rate is zero at the ghost layer
            WFP_weights(ghost_layer_cell_id) = 0.0;
            // the well rate is a sum of rates of ghost and top collector layers
            WFP_weights(top_collector_cell_id) = RFP_weights(ghost_layer_cell_id) + RFP_weights(top_collector_cell_id);
        }

        // void set_P_top(RealType rate)
        // {
        //     RealType factor{TwoPi / fluid.viscosity};
        //     RealType P_top = (rate / factor - (permeability * cell_volumes * (fluid.density * Gravity::value() * mesh_nodes - ext_pressure) * is_permeable.log_vals).sum() / std::log(R_ext / r_col)) /
        //                      ((permeability * cell_volumes * is_permeable.log_vals).sum() / std::log(R_ext / r_col));
        // }

        StepPropertyContainer get_RFP(
            RealType rate,
            RealType pressure = std::numeric_limits<double>::quiet_NaN()) const override
        {
            if (std::isnan(rate))
            { // define rate from pressure
                assert(!std::isnan(pressure));
                rate = 2 * std::numbers::pi / fluid.viscosity / log_dist * pressure * weights_sum;
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                pressure = rate/(2 * std::numbers::pi / fluid.viscosity / log_dist * weights_sum);
            }
            else
                assert("Incorrect injector regime!");

            assert(pressure > 0.0);
            assert(rate > 0.0);
            return ((rate / weights_sum) * RFP_weights).eval();
            //  {
            //     Logs::StepPropertyGrid{
            //         Logs::StepProperty{(RFP_weights * (rate / RFP_weights.sum())).eval()},
            //         grid},
            //     is_permeable};
        }

        StepPropertyContainer get_WFP(
            RealType rate,
            RealType pressure = std::numeric_limits<double>::quiet_NaN()) const
        {
            if (std::isnan(rate))
            { // define rate from pressure
                assert(!std::isnan(pressure));
                rate = 2 * std::numbers::pi / fluid.viscosity / log_dist * pressure * weights_sum;
            }
            else if (std::isnan(pressure))
            { // define pressure from rate
                assert(!std::isnan(rate));
                pressure = rate/(2 * std::numbers::pi / fluid.viscosity / log_dist * weights_sum);
            }
            else
                assert("Incorrect injector regime!");

            assert(pressure > 0.0);
            assert(rate > 0.0);

            return ((rate / weights_sum) * WFP_weights).eval();
        }

        const ptrdiff_t top_collector_cell_id{-1ll};
        const ptrdiff_t ghost_layer_cell_id{-1ll};

    protected:
        const StepPropertyContainer RFP_weights;
        StepPropertyContainer WFP_weights;
        RealType weights_sum;
        const PhaseProperties fluid;

    private:
        static ptrdiff_t layer_id(const auto &indicator)
        {
            const auto perforated_it = std::ranges::find(indicator.log_vals, 1.0);
            return std::distance(indicator.log_vals.cbegin(), perforated_it);
        }
        const RealType log_dist;
    };
} // GPN
