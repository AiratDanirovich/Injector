#pragma once

#include <memory>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/ConcreteGrids.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/Factory.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/InitialCondition.hpp>

#include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Properties/FlowField.hpp>

#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well.hpp>

namespace GPN
{
    namespace EqSolver
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

        template<typename Grid2D_t>
        struct SolverFactory
            : 
              public Logs::HydrodynamicLogsFactory
        {
            SolverFactory(RealType val, const Box &box, ptrdiff_t n1, ptrdiff_t n2)
                : SolverFactory{val, Grids::CylinderGridFactory::create(box, n1, n2)}
            {
            }

            /// @brief
            /// @param val initial condition const-value
            /// @param grid_factory
            SolverFactory(RealType val, const Grid2D_t &grid2D)
                : grid2D{grid2D},
                  Logs::HydrodynamicLogsFactory{grid2D->first_coord},
                  val{val},
                  well{water(), is_permeable, permeability}
            {
            }

            Properties::HeatConductivity conductivity_field() const
            {
                return {conductivity_field(Logs::RawdataFactory::generate_conductivity_StepProperty(
                    grid2D->first_coord.dual_stencils))};
            }

            Properties::HeatConductivity conductivity_field(const std::vector<RealType> &vals) const
            {
                return {
                    Logs::HeatConductivity{
                        Logs::StepPropertyGrid{
                            Logs::StepProperty{
                                vals},
                            grid2D->first_coord}},
                    grid2D};
            }

            Properties::HeatVolumetricCapacity capacity_field() const
            {
                auto solid_density_stencils{
                    Logs::RawdataFactory::generate_solid_density_StepProperty(
                        z_grid_stencils())};
                auto solid_density{
                    Logs::SolidDensity{
                        Logs::StepPropertyGrid{
                            Logs::StepProperty{
                                solid_density_stencils},
                            z_grid()}}};

                auto solid_specific_heatcapacity_stencils{
                    Logs::RawdataFactory::generate_solid_specific_heatcapacity_StepProperty(
                        z_grid_stencils())};
                auto solid_specific_heatcapacity{
                    Logs::SolidSpecificHeatCapacity{
                        Logs::StepPropertyGrid{
                            Logs::StepProperty{
                                solid_specific_heatcapacity_stencils},
                            z_grid()}}};

                auto solid_volumetric_heatcapacity{
                    Logs::SolidVolumetricHeatCapacity{
                        solid_density, solid_specific_heatcapacity}};

                auto capacity{
                    Logs::HeatVolumetricCapacity{
                        porosity, solid_volumetric_heatcapacity, water()}};
                return {
                    capacity,
                    grid2D};
            }

            Water water() const
            {
                return Phases::FluidFactory::create_water(300, 10);
            }

            const Well_KH &well_kh() const
            {
                return well;
            }

            Properties::ReservoirFlowField flow_field(RealType well_rate) const
            {
                return {
                    well_rate, well_kh(), *grid2D};
            }

            Problem::InitialCondition initial_state() const
            {
                return {
                    State::State2D::FillWithConst(
                        *grid2D, val)};
            }

            BoundaryConditions::BoundaryConditions boundary_conditions()
            {
                return {
                    *grid2D,
                    std::make_shared<BCFunctor>(BCFunctor{val})};
            }

            const Grids::GridDualStencils &z_grid_stencils() const
            {
                return grid2D->first_coord.dual_stencils;
            }
            const Grids::AxesGrid<CoordinateTypes::Z> &z_grid() const
            {
                return grid2D->first_coord;
            }

            Well_KH well;

            RealType val;

            const Grid2D_t grid2D;
        };

    } // EqSolver

} // GPN