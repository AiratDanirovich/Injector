#pragma once

namespace GPN
{
    struct InjectorRegimes
    {
        enum Type
        {
            FixedRate,
            FixedBottomHolePressure
        };
    };
} // GPN