#pragma once

#include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/History/InjectorRegimes.hpp>

#include <Injector/Properties/FlowField.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <
            typename Grid2D_t, typename Well_t,
            typename History_t, typename Fluid_t,
            typename Hydrodynamics_t>
        struct CompressibleRatesFactory
        {
            using OriginalGrid = typename Grid2D_t::OriginalGrid;

            CompressibleRatesFactory(
                cptr<Hydrodynamics_t> pressure_field,
                const cptr<Grid2D_t> grid2D_rocks,
                const Well_t &well,
                const ptr<History_t> history,
                const Fluid_t &fluid)
                : grid2D_rocks{grid2D_rocks},
                  grid2D{*(grid2D_rocks->grid2D)}, 
                  well{well},
                  history{history},
                  fluid{fluid},
                  pressure_field{pressure_field},
                  first_size{
                    grid2D_rocks->grid2D->first_coord().mesh_size()}, 
                  second_size{
                    grid2D_rocks->grid2D->second_coord().mesh_size()}
            {
            }

            /// @brief set the flow field at the next time moment
            /// @param t current time moment
            /// @param t_step step to the next time moment
            void set_flow_field(
                double t, RealType t_step)
            {
                // this method only works at FixedRate injection
                assert(history->regime() == InjectorRegimes::FixedRate);

                // the filed is updated at every time step
                // non-stationary hydrodynamics is assumed
                pressure_field->set_pressure_field(t_step, get_history_record());

                // volumetric flow field in two directions is calculated,
                // once the pressure field is calculated
                heat_flow_field =
                    std::make_shared<FaceProperties::HeatFlowField>(
                        // create ReservoirFlowField
                        FaceProperties::FlowFactory::create_from_well(
                            history->get_current_record(), well, *grid2D_rocks),
                        // multiple by heat capaity
                        fluid.volumetric_heat_capacity);

                const FaceValuesContainer axes1_value{
                    FaceValuesContainer::Zero(
                        grid2D->first_coord().dual_size(),
                        grid2D->second_coord().mesh_size())};
                FaceValuesContainer axes2_value{
                    FaceValuesContainer::Zero(
                        grid2D->first_coord().mesh_size(),
                        grid2D->second_coord().dual_size())};

                axes2_value.leftCols(Grid2D_t::l_margin+1ll).colwise() = well.get_RFP(get_history_record());
                
            //    axes2_value.middleCols(second_size - Grid2D_t::l_margin) = rock_P;

                heat_flow_field =
                    std::make_shared<FaceProperties::HeatFlowField>(
                        ReservoirFlowField{
                            axes1_value,
                            axes2_value},
                        // multiple by heat capaity
                        fluid.volumetric_heat_capacity);
            }

        public:
            const OriginalGrid &grid2D;
            const cptr<Grid2D_t> grid2D_rocks;
            const Well_t &well;
            const ptr<History_t> history;
            const Fluid_t &fluid;
            ptr<Hydrodynamics_t> pressure_field;

        protected:
            const ptrdiff_t first_size, second_size;

            cptr<FaceProperties::HeatFlowField> heat_flow_field;
            //    cptr<FaceProperties::ReservoirFlowField> volumetric_flow_field;

            const auto get_history_record() const
            {
                return history->get_current_record();
            }
        };

    } // FaceProperties
} // GPN
