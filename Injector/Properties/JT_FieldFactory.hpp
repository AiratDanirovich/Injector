#pragma once

#include <cassert>
#include <utility>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/PhysicalField.hpp>

namespace GPN
{
    namespace Properties
    {
        struct JT_FieldFactory
        {
            static auto create_spatial(
                const auto &rates_factory)
            {
                const auto /*not &*/ grid2D_ptr{rates_factory.grid2D};
                const auto first_coord_size{
                    grid2D_ptr->first_coord().mesh_size()};
                const auto second_coord_size{
                    grid2D_ptr->second_coord().mesh_size()};
                GridNodeValues2D values{
                    GridNodeValues2D::Zero(
                        first_coord_size,
                        second_coord_size)};

                const auto &flux2_pos{rates_factory.get_heat_flow_in_axes2_pos()};
                const auto &flux2_neg{rates_factory.get_heat_flow_in_axes2_neg()};
                const auto &permeability{rates_factory.pressure_field->permeability};
                const auto &pressure{rates_factory.get_pressure_field().values()};

                assert(flux2_pos.cols() == flux2_neg.cols());
                assert(flux2_pos.cols() == second_coord_size + 1ll);

                // row-array
                MeshNodesContainerT buffer{MeshNodesContainer::Zero(second_coord_size)};
                for (auto row{0ll}; row < first_coord_size; ++row)
                {
                    if (permeability(row) > 0.0)
                    {
                        buffer =
                            (flux2_pos.row(row).head(second_coord_size) -
                             flux2_neg.row(row).tail(second_coord_size)) *
                            pressure.row(row);

                        buffer.head(second_coord_size - 1ll) +=
                            flux2_neg.block(row, 1ll, 1ll, second_coord_size - 1ll) *
                            pressure.row(row).tail(second_coord_size - 1ll);
                        buffer.tail(second_coord_size - 1ll) -=
                            flux2_pos.block(row, 1ll, 1ll, second_coord_size - 1ll) *
                            pressure.row(row).head(second_coord_size - 1ll);

                        values.row(row) = buffer*rates_factory.fluid.JT;
                    }
                }
                
                for(auto row{0ll}; row < values.rows(); ++row)
                    for(auto col{0ll}; col < values.cols(); ++col)
                        assert(!std::isnan(values(row,col)) && !std::isinf(values(row,col)));

                return JT_SpatialComponent{std::move(values), grid2D_ptr};
            }
            
            static auto create_temporal(
                const auto &rates_factory)
            {
                const auto& p_field{*(rates_factory.pressure_field)};
                const auto& porosity{p_field.porosity.log_vals};
                const auto& adiabatic_factor{
                    p_field.fluid.adiabatic_factor};

                const auto /*not &*/ grid2D_ptr{rates_factory.grid2D};
                const auto first_coord_size{
                    grid2D_ptr->first_coord().mesh_size()};
                const auto second_coord_size{
                    grid2D_ptr->second_coord().mesh_size()};
                GridNodeValues2D values{
                    p_field.delta_pressure()*(grid2D_ptr->volumes().colwise()*(porosity*adiabatic_factor))};

                // temporal contribution of JT
                // inside the sandface is assumed zero
                values.leftCols(rates_factory.grid2D_rocks->l_margin) = 0.0;

                for(auto row{0ll}; row < values.rows(); ++row)
                    for(auto col{0ll}; col < values.cols(); ++col)
                        assert(!std::isnan(values(row,col)) && !std::isinf(values(row,col)));

                return JT_SpatialComponent{std::move(values), grid2D_ptr};
            }
        };
    } // Properties
} // GPN