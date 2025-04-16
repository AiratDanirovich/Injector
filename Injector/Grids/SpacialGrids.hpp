#pragma once

#include <Injector/Grids/Grids1D.hpp>


namespace GPN
{
    namespace Grids
    {
        /// @brief Stencils of spacial grid.
        /// It consists of dual nodes -- nodes of property jump change
        struct SpacialGridStencils : public GridDualStencils
        {
            using GridDualStencils::GridDualStencils;
        };

        /// @brief Refined spacial grid
        struct SpacialGrid : public GridDual
        {
            using GridDual::GridDual;
        };

    } // Grids
} // GPN
