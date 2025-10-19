#pragma once

#include <memory>

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/BoundaryConditions.hpp>

namespace GPN
{
    namespace Heat
    {
        template <typename History_t>
        struct HeatBC : public BoundaryConditions::GeneralBC
        {
            template <typename Grid2D_t>
            HeatBC(const cptr<Grid2D_t> &grid,
                   cptr<const BCFunctorBase> functor,
                   const cptr<History_t> history)
                : BoundaryConditions::GeneralBC{grid, functor, BoundaryCondition::BCType::second},
                history{history}
            {
            }

            void set_bc_type(const RealType t)
            {
                this->t = t;
            //    bc_types[south_id] = BoundaryCondition::BCType::second; // top boundary
            //    bc_types[north_id] = BoundaryCondition::BCType::second; // bottom boundary
            //    bc_types[west_id] = BoundaryCondition::BCType::second;  // well axis of symmetry
                bc_types[east_id] = BoundaryCondition::BCType::second;  // external contour
            }

        protected:
            const cptr<History_t> history;
        };
    } // Heat
} // GPN