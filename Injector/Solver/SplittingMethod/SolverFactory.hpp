#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Grids/PhysicalField.hpp>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/Factory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <Injector/Solver/SplittingMethod/Solver.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace SplittingMethod
        {
            struct Factory
            {
                struct BCFunctor : public GPN::BoundaryConditions::BCFunctorBase
                {
                    BCFunctor(RealType val) : val{val} {}

                    RealType operator()(
                        RealType x, RealType y, RealType t) const override
                    {
                        return val;
                    }

                protected:
                    RealType val;
                };

                static auto make_solver(RealType val)
                {
                    using namespace GPN;
                    using namespace GPN::Grids;
                    using namespace GPN::Logs;
                    using namespace GPN::EqSolver;
                    using namespace GPN::EqSolver::State;
                    using namespace GPN::EqSolver::Problem;
                    using namespace GPN::EqSolver::SplittingMethod;

#pragma region GRID_2D
                    const cptr<StructuredCylinderGrid2DAxisymmetric> grid2D{
                        std::make_shared<StructuredCylinderGrid2DAxisymmetric>(
                            Grids::Factory::create_cylinder_grid_2D(
                                Box{Segment{0, 1}, Segment{0, 1}}, 301, 501))};
                    const auto &z_grid_stencils{grid2D->first_coord.dual_stencils};
                    const auto &z_grid{grid2D->first_coord};
#pragma endregion
#pragma region HEAT-CONDUCTIVITY
                    // generate heat conductivity field
                    Properties::HeatConductivity conductivity_field{
                        Properties::Factory::generate_heatconductivity_Property(
                            Logs::Factory::generate_conductivity_StepProperty(
                                grid2D->first_coord.dual_stencils),
                            grid2D)};
#pragma endregion
#pragma region WATER
                    auto water{Phases::FluidFactory::create_water(300, 10)};
#pragma endregion
#pragma region HEAT-VOLUMETRIC-CAPACITY
                    Properties::HeatVolumetricCapacity capacity_field{
                        Properties::Factory::generate_volumetric_heatcapacity_Property(
                            Logs::Factory::generate_is_permeable_StepProperty(z_grid_stencils),
                            Logs::Factory::generate_porosity_StepProperty(z_grid_stencils),
                            Logs::Factory::generate_solid_density_StepProperty(
                                z_grid_stencils),
                            Logs::Factory::generate_solid_specific_heatcapacity_StepProperty(
                                z_grid_stencils),
                            water,
                            grid2D)};
#pragma endregion

                    InitialCondition
                        initial_state{
                            State2D::FillWithConst(
                                *grid2D, val)};

                    GPN::BoundaryConditions::BoundaryConditions bc{
                        *grid2D, std::make_shared<BCFunctor>(BCFunctor{val})};

                    const double tol = 1E-8;

                    return Solver{
                        conductivity_field, grid2D,
                        capacity_field,
                        initial_state,
                        bc, 0.0};
                }
            };
        } // SplittingMethod
    } // EqSolver
} // GPN
