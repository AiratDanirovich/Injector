#pragma once

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/BoundaryConditions.hpp>

namespace GPN
{
    namespace Heat
    {
        struct HeatBC : public BoundaryConditions::GeneralBC
        {
            using BoundaryConditions::GeneralBC::GeneralBC;

            void set_bc_type(const RealType t){
                this->t = t;
                bc_types[south_id] = BoundaryCondition::BCType::second; // top boundary
                bc_types[north_id] = BoundaryCondition::BCType::second; // bottom boundary
                bc_types[west_id] = BoundaryCondition::BCType::second; // well axis of symmetry
                bc_types[east_id] = BoundaryCondition::BCType::second; // external contour
            }

        };
    } // Heat
} // GPN