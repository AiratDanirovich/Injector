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

#include <Injector/Model/Hydrodynamic/SomeRatesFactory.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <
            typename Grid2D_t, typename Well_t,
            typename History_t, typename Fluid_t,
            typename Hydrodynamics_t>
        struct IncompressibleRatesFactory : public SomeRatesFactory<Grid2D_t, Well_t, History_t, Fluid_t, Hydrodynamics_t>
        {
            using Base = SomeRatesFactory<Grid2D_t, Well_t, History_t, Fluid_t, Hydrodynamics_t>;
            using Base::get_history_record;
            using Base::history;
            using Base::solution;

            IncompressibleRatesFactory(
                cptr<Hydrodynamics_t> pressure_field,
                const cptr<Grid2D_t> grid2D_rocks,
                const Well_t &well,
                const ptr<History_t> history,
                const Fluid_t &fluid)
                : Base{pressure_field, grid2D_rocks, well, history, fluid},
                  fluid{fluid},
            {
            }

            /// @brief set the flow field at the next time moment
            /// @param t current time moment
            /// @param t_step step to the next time moment
            void set_flow_field(
                double t, RealType t_step)
            {
                Base::set_flow_field(
                    t, t_step,
                    [&well->rfp](const ptrdiff_t /*count*/, const ptrdiff_t /*col*/)
                    {
                        return well->rfp
                    });
            }

            const auto get_spatial_JT_contribution() const
            {
                return Properties::JT_FieldFactory::create_spatial(*this).values();
            }

            const auto get_temporal_JT_contribution() const
            {
                throw std::runtime_error("The temporal JT contribution was not tested for Incompressible fluid.");
                //    return Properties::JT_FieldFactory::create_temporal(*this).values();
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

        protected:
            cptr<FaceProperties::HeatFlowField> heat_flow_field;
            cptr<FaceProperties::ReservoirFlowField> volumetric_flow_field;
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