#pragma once
#include <memory>
#include <limits>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Properties/PhysicalField.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <typename Grid2D_t>
        struct ZeroRatesFactory
        {
            ZeroRatesFactory(
                const cptr<Grid2D_t> grid2D,
                const Logs::IsPermeable &is_permeable,
                const Logs::ExternalPressure &ext_pressure)
                : heat_flow_field{
                      std::make_shared<FaceProperties::HeatFlowField>(
                          FaceProperties::FlowFactory::zero_flow(
                              is_permeable, *grid2D),
                          1.0)},
                  P_ext{Properties::FieldFactory::create(ext_pressure, grid2D)},
                  grid2D{grid2D}
            {
            }

            void set_flow_field(double, RealType)
            {
            }
            
            const auto get_spatial_JT_contribution() const
            {
                return GridNodeValues2D::Zero(
                    grid2D->first_coord().mesh_size(),
                    grid2D->second_coord().mesh_size());
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

            const auto get_pressure_field() const
            {
                return P_ext;
            }

        public:
            const cptr<Grid2D_t> grid2D;
        protected:
            const cptr<FaceProperties::HeatFlowField> heat_flow_field;

            const Properties::Pressure<Grid2D_t> P_ext;
        };
    } // FaceProperties
} // GPN