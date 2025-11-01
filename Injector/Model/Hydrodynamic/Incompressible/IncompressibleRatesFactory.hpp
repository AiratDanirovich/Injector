#pragma once
#include <algorithm>
#include <cassert>
#include <memory>
#include <limits>
#include <vector>

#include <Injector/Grids/Defines.h>

#include <Injector/History/InjectorRegimes.hpp>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Properties/JT_FieldFactory.hpp>

#include <Injector/Solver/State2D.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <
            typename Grid2D_t, typename Well_t,
            typename History_t, typename Fluid_t,
            typename Hydrodynamics_t>
        struct IncompressibleRatesFactory
        {
            using OriginalGrid = typename Grid2D_t::OriginalGrid;

            struct Solution
            {
                Solution(const auto history) 
                {
                    times.reserve(history->size());
                    states.reserve(history->size());
                }
                std::vector<RealType> times;
                std::vector<EqSolver::State::State2D> states;
            };

            Solution solution;

            IncompressibleRatesFactory(
                cptr<Hydrodynamics_t> pressure_field,
                const cptr<Grid2D_t> grid2D_rocks,
                const Well_t &well,
                const ptr<History_t> history,
                const Fluid_t &fluid)
                : grid2D_rocks{grid2D_rocks},
                  grid2D{grid2D_rocks->grid2D},
                  well{well},
                  history{history},
                  fluid{fluid},
                  pressure_field{pressure_field},
                  solution{history}
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
                            history->get_current_record(), well, *grid2D),
                        // multiply by heat capacity
                        fluid.volumetric_heat_capacity);
            }

            const auto get_spatial_JT_contribution() const
            {
                return Properties::JT_FieldFactory::create(*this).values();
            }

            const auto get_heat_flow_in_axes1() const
            {
                return heat_flow_field->axes1_as_face_normal_pos + heat_flow_field->axes1_as_face_normal_neg;
            }
            const auto get_heat_flow_in_axes2() const
            {
                return heat_flow_field->axes2_as_face_normal_pos + heat_flow_field->axes2_as_face_normal_neg;
            }

            const auto &get_heat_flow_in_axes1_pos() const
            {
                return heat_flow_field->axes1_as_face_normal_pos;
            }
            const auto &get_heat_flow_in_axes2_pos() const
            {
                return heat_flow_field->axes2_as_face_normal_pos;
            }
            const auto &get_heat_flow_in_axes1_neg() const
            {
                return heat_flow_field->axes1_as_face_normal_neg;
            }
            const auto &get_heat_flow_in_axes2_neg() const
            {
                return heat_flow_field->axes2_as_face_normal_neg;
            }
            const auto &get_pressure_field() const
            {
                return pressure_field->current_pressure();
            }

        public:
            const cptr<OriginalGrid> grid2D;
            const cptr<Grid2D_t> grid2D_rocks;
            const Well_t &well;
            const ptr<History_t> history;
            const Fluid_t &fluid;
            ptr<Hydrodynamics_t> pressure_field;

        protected:
            cptr<FaceProperties::HeatFlowField> heat_flow_field;
            //    cptr<FaceProperties::ReservoirFlowField> volumetric_flow_field;
            
            const auto get_history_record() const
            {
                return history->get_current_record();
            }
        };

        template <typename Grid2D_t, typename Fluid_t>
        struct HorizontalRatesFactory
        {
            HorizontalRatesFactory(
                const RealType well_rate,
                const Fluid_t &fluid,
                const cptr<Grid2D_t> grid2D,
                const Logs::IsPermeable &is_permeable)
                : heat_flow_field{
                      std::make_shared<FaceProperties::HeatFlowField>(
                          FaceProperties::FlowFactory::horizontal_flow(
                              well_rate, is_permeable, *grid2D),
                          fluid.volumetric_heat_capacity)}
            {
            }

            void set_flow_field(
                double, RealType)
            {
            }

            const auto &get_heat_flow_in_axes1() const
            {
                return heat_flow_field->axes1_as_face_normal;
            }
            const auto &get_heat_flow_in_axes2() const
            {
                return heat_flow_field->axes2_as_face_normal;
            }
            const auto &get_heat_flow_in_axes1_pos() const
            {
                return heat_flow_field->axes1_as_face_normal_pos;
            }
            const auto &get_heat_flow_in_axes2_pos() const
            {
                return heat_flow_field->axes2_as_face_normal_pos;
            }
            const auto &get_heat_flow_in_axes1_neg() const
            {
                return heat_flow_field->axes1_as_face_normal_neg;
            }
            const auto &get_heat_flow_in_axes2_neg() const
            {
                return heat_flow_field->axes2_as_face_normal_neg;
            }

        protected:
            const cptr<Grid2D_t> grid2D;
            const cptr<FaceProperties::HeatFlowField> heat_flow_field;
        };
    } // FaceProperties
} // GPN