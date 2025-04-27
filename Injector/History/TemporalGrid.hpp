#pragma once

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace CoordinateTypes
    {
        struct Time : public GeneralCoordinate
        {
        };
    } // CoordinateTypes

    namespace Grids
    {
        /// @brief Stencils of temporal grid consists of dual nodes -- moments of rate change
        struct TemporalGridDualStencils : public GridDualStencils
        {
            using GridDualStencils::GridDualStencils;
        };

        /// @brief Refined temporal grid
        struct TemporalGridDual : public GridDual
        {
            using GridDual::GridDual;
        };

    } // Grids
} // GPN
