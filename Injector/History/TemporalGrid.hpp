#pragma once

#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace CoordinateTypes
    {
        struct Time : public GeneralCoordinate
        {
            /// @brief 
            /// @param nodes Nodes of dual mesh
            /// @return Centers of control volumes
            static auto cell_centers(const DualNodesContainer& nodes)
            {
                assert(nodes.size() > 2ll);

                const RealType tol = 1e-12;

                auto size{nodes.size()-1ll};
                MeshNodesContainer out(size);
                // centers of boundary cells are moved to the domain boundary
                for(auto idx{0ll}; idx < size; ++idx)
                    out(idx) = (nodes(idx+1) + nodes(idx))/2.0;

                return out;
            }

        };
    } // CoordinateTypes

    namespace Grids
    {
        /// @brief Stencils of temporal grid
        /// consists of dual nodes -- moments of rate change
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
