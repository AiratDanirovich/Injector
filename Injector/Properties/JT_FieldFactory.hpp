#pragma once

#include <cassert>

#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/History/RatesFactory.hpp>

namespace GPN
{
    namespace Properties
    {
        struct JT_FieldFactory
        {
            static auto create(
                const auto &rates_factory)
            {
                const auto /*&*/ grid2D_ptr{rates_factory.grid2D};
                const auto first_coord_size{
                    grid2D_ptr->first_coord.mesh_size()};
                const auto second_coord_size{
                    grid2D_ptr->second_coord.mesh_size()};
                GridNodeValues2D values{
                    GridNodeValues2D::Zero(
                        first_coord_size,
                        second_coord_size)};

                const auto &flux2_pos{rates_factory.get_heat_flow_in_axes2_pos()};
                const auto &flux2_neg{rates_factory.get_heat_flow_in_axes2_neg()};
                const auto &permeability{rates_factory.pressure_field->permeability};
                const auto &pressure{rates_factory.get_pressure_field().its_values};

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

                        values.row(row) = buffer;
                    }
                }

                return JT_SpatialComponent{values, grid2D_ptr};
            }
        };
    } // Properties
} // GPN