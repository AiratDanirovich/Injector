#pragma once
#include <algorithm>
#include <cassert>
#include <memory>
#include <limits>

#include <Injector/History/History.hpp>
#include <Injector/Properties/FlowField.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Properties/JT_FieldFactory.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <
            typename Grid2D_t,
            typename Well_t,
            typename Fluid_t,
            typename Hydrodynamics_t>
        struct IncompressibleRatesFactory
        {
            IncompressibleRatesFactory(
                cptr<Hydrodynamics_t> pressure_field,
                const cptr<Grid2D_t> grid2D,
                const Well_t &well,
                const History &history,
                const Fluid_t &fluid)
                : grid2D{grid2D},
                  well{well},
                  history{history},
                  fluid{fluid},
                  pressure_field{pressure_field},
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

                    // volumetric flow field in two directions
                    heat_flow_field =
                        std::make_shared<FaceProperties::HeatFlowField>(
                            // create ReservoirFlowField
                            FaceProperties::FlowFactory::create_from_well(
                                history.get_record(pos), well, *grid2D),
                            // multiple by heat capaity
                            fluid.volumetric_heat_capacity);
                }

                pressure_field->set_pressure_field(t_step, get_history_record());
            }

            const auto get_spatial_JT_contribution() const
            {
                return Properties::JT_FieldFactory::create(*this).values();
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

            const auto& get_pressure_field() const
            {
                return pressure_field->current_pressure();
            }

        public:
            const cptr<Grid2D_t> grid2D;
            const Well_t &well;
            const History &history;
            const Fluid_t &fluid;
            ptr<Hydrodynamics_t> pressure_field;

        protected:
            ptr<FaceProperties::HeatFlowField> heat_flow_field;

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