#pragma once

#include <Injector/Model/Completion.hpp>
#include <Injector/Model/Well/Well.hpp>

namespace GPN
{
    struct WellHolesFactory
    {
        static WellHoles create(const Completion::Casing<Completion::VarRing>& completion)
        {
            return WellHoles{
                TubeInnerRadius{completion.flow_radius(0ll)}, 
                ColumnOuterRadius{completion.column_outer_radius(0ll)},
                SandfaceRadius{completion.sandface_radius(0ll)}};
        }
    };

} // GPN