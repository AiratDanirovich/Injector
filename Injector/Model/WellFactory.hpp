#pragma once

#include <Injector/Model/Completion.hpp>
#include <Injector/Model/Well.hpp>

namespace GPN
{
    struct WellHolesFactory
    {
        static WellHoles create(const Completion::Casing& completion)
        {
            return WellHoles{
                TubeInnerRadius{completion.flow_radius}, 
                ColumnOuterRadius{completion.column_outer_radius},
                SandfaceRadius{completion.sandface_radius}};
        }
    };

} // GPN