#pragma once

#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/History/RatesFactory.hpp>

namespace GPN
{
    namespace Properties
    {
        struct JT_FieldFactory
        {
            static auto create(
                const auto& flow)
            {
                const auto grid2D_ptr{flow.grid2D};
                GridNodeValues2D values;


                return JT_SpatialComponent{values, grid2D_ptr};
            }
        };
    } // Properties
} // GPN