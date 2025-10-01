#pragma once
#include <algorithm>
#include <cassert>
#include <memory>
#include <limits>

#include <Injector/History/History.hpp>
#include <Injector/Properties/FlowField.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <typename Grid2D_t, typename Well_t, typename Fluid_t>
        struct RatesFactory
        {
            RatesFactory(
                const cptr<Grid2D_t> grid2D,
                const Well_t &well,
                const History &history,
                const Fluid_t &fluid)
                : grid2D{grid2D},
                  well{well},
                  history{history},
                  fluid{fluid},
                  pos{-1ll}
            {
            }

            void set_flow_field(
                double t, RealType t_step)
            {
                // first history.time_moment greater than t_mid
                const auto it{std::upper_bound(
                    history.time_moments.cbegin(),
                    history.time_moments.cend(),
                    t + t_step / 2.0)};
                // corresponding position in history.rates
                const auto pos_new{std::distance(history.time_moments.cbegin(), it) - 1ll};
                if (pos_new > pos)
                { // the filed is only updated if a new history interval is set
                    pos = pos_new;
                    assert(pos >= 0ll);
                    // this method only works at FixedRate injection
                    assert(history.regimes[pos] == InjectorRegimes::FixedRate);

                    // volumetric flow filed in two directions
                    volumetric_flow_field = 
                        std::make_shared<FaceProperties::ReservoirFlowField>(
                            FaceProperties::FlowFactory::create_from_well(
                                history.get_record(pos), well, *grid2D));
                    // heat flow field in two directions.
                    // First, the flow field is copied
                    heat_flow_field =
                        std::make_shared<FaceProperties::ReservoirFlowField>(*volumetric_flow_field);
                    // Second, the flow files is multiplied by volumetric heat capacity of fluid
                    FaceProperties::multiply(*heat_flow_field, fluid.volumetric_heat_capacity);
                }
            }

            const auto &get_flow_in_axes1() const
            {
                return heat_flow_field->axes1_as_face_normal;
            }
            const auto &get_flow_in_axes2() const
            {
                return heat_flow_field->axes2_as_face_normal;
            }

            const auto get_temperature() const
            {
                return history.temps(pos);
            }
            const auto get_rate() const
            {
                return history.rates(pos);
            }
            const auto get_pressure() const
            {
                return history.pressure(pos);
            }
            const auto get_history_record() const
            {
                return history.get_record(pos);
            }

        protected:
            const cptr<Grid2D_t> grid2D;
            const Well_t &well;
            const History &history;
            const Fluid_t &fluid;
            cptr<FaceProperties::ReservoirFlowField> heat_flow_field;
            cptr<FaceProperties::ReservoirFlowField> volumetric_flow_field;

        private:
            std::ptrdiff_t pos{-1ll};
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
                      std::make_shared<FaceProperties::ReservoirFlowField>(
                          FaceProperties::FlowFactory::horizontal_flow(
                              well_rate, is_permeable, *grid2D))}
            {
                FaceProperties::multiply(*heat_flow_field, fluid.volumetric_heat_capacity);
            }

            void set_flow_field(
                double, RealType)
            {
            }

            const auto &get_flow_in_axes1() const
            {
                return heat_flow_field->axes1_as_face_normal;
            }
            const auto &get_flow_in_axes2() const
            {
                return heat_flow_field->axes2_as_face_normal;
            }

        protected:
            const cptr<Grid2D_t> grid2D;
            const cptr<FaceProperties::ReservoirFlowField> heat_flow_field;
        };

        template <typename Grid2D_t>
        struct ZeroRatesFactory
        {
            ZeroRatesFactory(
                const cptr<Grid2D_t> grid2D,
                const Logs::IsPermeable &is_permeable)
                : heat_flow_field{
                      std::make_shared<FaceProperties::ReservoirFlowField>(
                          FaceProperties::FlowFactory::zero_flow(
                              is_permeable, *grid2D))}
            {
            }

            void set_flow_field(double, RealType)
            {
            }

            const auto &get_flow_in_axes1() const
            {
                return heat_flow_field->axes1_as_face_normal;
            }
            const auto &get_flow_in_axes2() const
            {
                return heat_flow_field->axes2_as_face_normal;
            }

        protected:
            const cptr<FaceProperties::ReservoirFlowField> heat_flow_field;
        };
    } // Properties
} // GPN