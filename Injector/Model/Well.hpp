#pragma once

#include <numbers>

#include <Eigen/Core>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/History/History.hpp>

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

    struct Well
    {
        template <typename IsPermeable_t, typename Permeability_t, typename ExternalPressure_t>
        Well(
            RealType R_ext,
            RealType r_col,
            RealType r_tube,
            const PhaseProperties &fluid,
            const Friction &friction,
            const IsPermeable_t &is_permeable,
            const Permeability_t &permeability,
            const ExternalPressure_t &ext_pressure)
            : R_ext{R_ext},
              r_col{r_col},
              r_tube{r_tube},
              fluid{fluid},
              friction{friction},
              is_permeable{is_permeable.log_vals},
              permeability{permeability.log_vals},
              ext_pressure{ext_pressure.log_vals},
              mesh_nodes{is_permeable.grid.get_mesh_nodes()},
              cell_volumes{is_permeable.grid.get_dual_steps()}
        {
        }

        void set_P_top(RealType rate)
        {
            RealType factor{TwoPi / fluid.viscosity};
            RealType P_top = (rate / factor - (permeability * cell_volumes * (fluid.density * Gravity::value() * mesh_nodes - ext_pressure) * is_permeable).sum() / std::log(R_ext / r_col)) /
                             ((permeability * cell_volumes * is_permeable).sum() / std::log(R_ext / r_col));

                             
        }

    protected:
        const RealType R_ext, r_col, r_tube;
        const PhaseProperties fluid;
        const Logs::StepPropertyContainer
            is_permeable,
            permeability,
            ext_pressure,
            mesh_nodes,
            cell_volumes;
        const Friction friction;

    private:
        RealType TwoPi{2.0*std::numbers::pi};
        //     RealType P_top;
    };

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
