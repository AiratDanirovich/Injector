#pragma once

#include <numbers>
#include <cassert>

#include <Eigen/Core>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Properties/Logs.hpp>
// #include <Injector/History/History.hpp>

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
            const RealType well_radius)
            : tube_radius{tube_radius},
              well_radius{well_radius}
        {
        }

        const RealType tube_radius;
        const RealType well_radius;
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

    //         WellRadius well_radius,
    //         const Permeability& permeability) noexcept
    //     : well_radius{well_radius}
    //     , permeability{permeability}
    //     {}

    // protected:
    //     WellRadius well_radius;
    //     Permeability permeability;
    // };

} // GPN
