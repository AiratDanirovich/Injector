#pragma once

#include <vector>
#include <numbers>
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

    /// @brief Descriptor for the well circular desing,
    /// contains radii of tube < column < sandface.
    /// Cement is between column and sandface
    struct WellHoles
    {
        WellHoles(
            const RealType tube_radius,
            const RealType column_radius,
            const RealType sandface_radius)
            : tube_radius{tube_radius},
              column_radius{column_radius},
              sandface_radius{sandface_radius}
        {
            assert(tube_radius < sandface_radius);
        }

        std::vector<RealType> generate_uniform_radial_grid(
            const RealType r_min, const RealType r_max, const ptrdiff_t r_nodes)
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
        const RealType column_radius;
        const RealType sandface_radius;
    };

    /// @brief Descriptor of materials that fill the
    /// rings that form the well up to sandface.
    /// There is a variation along the verticle direction.
    /// Tube ends at the depth z_tube.
    struct WellMaterial
    {
        WellMaterial(
            const WellHoles &well_holes,
            const RealType tube_lambda,
            const RealType column_lambda,
            const RealType sandface_lambda,
            const RealType z_tube)
            : well_holes{well_holes},
              tube_lambda{tube_lambda},
              column_lambda{column_lambda},
              sandface_lambda{sandface_lambda},
              z_tube{z_tube}
        {
        }

        /// @brief Depth of the tube
        const RealType z_tube;

        /// @brief well concentric geometry
        const WellHoles well_holes;
        const RealType tube_lambda;
        const RealType column_lambda;
        const RealType sandface_lambda;
    };

    struct IWellDesign
    {
        IWellDesign(const Logs::IsPermeable &is_permeable)
            : is_permeable{is_permeable}
        {
        }
        using Grid_t = Logs::StepPropertyGrid::Grid_t;
        virtual LogValuesContainer get_RFP(RealType rate) const = 0;

        const Logs::IsPermeable is_permeable;
    };

    struct Well_KH : public IWellDesign
    {
        Well_KH(
            const PhaseProperties &fluid,
            const Logs::IsPermeable &is_permeable,
            const StepPropertyContainer &permeability)
            : IWellDesign{is_permeable},
              temp{permeability * is_permeable.grid.dual_steps * (StepPropertyContainer)is_permeable}
        {
            assert(permeability.size() == is_permeable.grid.dual_steps.size());
            assert(is_permeable.size() == is_permeable.grid.dual_steps.size());

            temp_sum = temp.sum();
        }

        // void set_P_top(RealType rate)
        // {
        //     RealType factor{TwoPi / fluid.viscosity};
        //     RealType P_top = (rate / factor - (permeability * cell_volumes * (fluid.density * Gravity::value() * mesh_nodes - ext_pressure) * is_permeable.log_vals).sum() / std::log(R_ext / r_col)) /
        //                      ((permeability * cell_volumes * is_permeable.log_vals).sum() / std::log(R_ext / r_col));
        // }

        StepPropertyContainer get_RFP(RealType rate) const override
        {
            return ((rate / temp_sum) * temp).eval();
            //  {
            //     Logs::StepPropertyGrid{
            //         Logs::StepProperty{(temp * (rate / temp.sum())).eval()},
            //         grid},
            //     is_permeable};
        }

    protected:
        //    const PhaseProperties fluid;
        const StepPropertyContainer temp;
        RealType temp_sum;
        //    const Logs::IsPermeable is_permeable;
    };

    // struct Well : public IWellDesign
    // {

    //     template <typename IsPermeable_t, typename Permeability_t, typename ExternalPressure_t>
    //     Well(
    //         RealType R_ext,
    //         RealType r_col,
    //         RealType r_tube,
    //         const PhaseProperties &fluid,
    //         const Friction &friction,
    //         const IsPermeable_t &is_permeable,
    //         const Permeability_t &permeability,
    //         const ExternalPressure_t &ext_pressure)
    //         : R_ext{R_ext},
    //           r_col{r_col},
    //           r_tube{r_tube},
    //           fluid{fluid},
    //           friction{friction},
    //           is_permeable{is_permeable},
    //           permeability{permeability.log_vals},
    //           ext_pressure{ext_pressure.log_vals},
    //           mesh_nodes{is_permeable.grid.get_mesh_nodes()},
    //           cell_volumes{is_permeable.grid.get_dual_steps()}
    //     //      ,
    //     //      grid{is_permeable.grid}
    //     {
    //     }

    //     void set_P_top(RealType rate)
    //     {
    //         RealType factor{TwoPi / fluid.viscosity};
    //         RealType P_top = (rate / factor - (permeability * cell_volumes * (fluid.density * Gravity::value() * mesh_nodes - ext_pressure) * is_permeable.log_vals).sum() / std::log(R_ext / r_col)) /
    //                          ((permeability * cell_volumes * is_permeable.log_vals).sum() / std::log(R_ext / r_col));
    //     }

    //     Logs::RFP get_RFP(RealType rate, const Grid_t &grid) const override
    //     {
    //         const auto temp{(permeability * cell_volumes * is_permeable.log_vals).eval()};

    //         return {
    //             Logs::StepPropertyGrid{
    //                 Logs::StepProperty{(temp * (rate / temp.sum())).eval()},
    //                 grid},
    //             is_permeable};
    //     }

    // protected:
    //     const RealType R_ext, r_col, r_tube;
    //     const PhaseProperties fluid;
    //     const Logs::StepPropertyContainer
    //         //    is_permeable,
    //         permeability,
    //         ext_pressure,
    //         mesh_nodes,
    //         cell_volumes;
    //     const Logs::IsPermeable is_permeable;
    //     //     const Grid_t &grid;
    //     const Friction friction;

    // private:
    //     RealType TwoPi{2.0 * std::numbers::pi};
    //     //     RealType P_top;
    // };

    // struct Well
    // {
    //     Well(

    //         WellRadius sandface_radius,
    //         const Permeability& permeability) noexcept
    //     : sandface_radius{sandface_radius}
    //     , permeability{permeability}
    //     {}

    // protected:
    //     WellRadius sandface_radius;
    //     Permeability permeability;
    // };

} // GPN
