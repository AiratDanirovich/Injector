#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>

namespace GPN
{
    namespace Model
    {
        struct Layers
        {
        };

        namespace Injector
        {

            struct FlowField
            {};


            // template <typename Phase_t = GPN::PhaseProperties>
            // struct Hydrodynamic2D
            // {

            //     Hydrodynamic2D(
            //         const Well &well,
            //         const Phase_t &fluid,
            //         const History &history,
            //         const Grid_t &grid) noexcept
            //         : history{history}
            //     {
            //     }

            //     template <Directions dir>
            //     auto flow_field() const
            //     {
            //         static_assert(
            //             dir < Dimensions::size,
            //             "Number of requsted coordinate direction exceeds "
            //             "the spacial problem dimension.");
            //         if constexpr (dir == Directions::x1)
            //         {
            //             FluxComponentContainer
            //                 out(
            //                     grid.first_coord.size(),
            //                     grid.second_coord.dual_size());
            //             out.setZero();
            //             return out;
            //         }
            //         if constexpr (dir == Directions::x2)
            //         {
            //             FluxComponentContainer
            //                 out(
            //                     grid.first_coord.dual_size(),
            //                     grid.second_coord.size());
            //             out.setZero();
            //             return out;
            //         }
            //     }

            // protected:
            //     Well well;
            //     const History &history;
            //     Phase_t fluid;
            //     const Grid_t &grid
            // };
        }
    }

} // GPN
