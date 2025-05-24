#pragma once

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Properties/FlowField.hpp>

namespace GPN
{
    namespace Problems
    {
        struct HeatProblem
        {
            BoundaryConditions::BoundaryCondition bc;
            EqSolver::Problem::InitialCondition ic;

            Properties::FlowField flow_field;

            Grids::StructuredCylinderGrid2DAxisymmetric grid;


        };

    } // Problems

} // GPN