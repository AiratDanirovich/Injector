#pragma once

#include <Injector/Model/Phases/PhaseProperties.hpp>

namespace GPN
{
    namespace Completion
    {
        struct MaterialType
        {
            enum
            {
                Flow = 0,
                Tube = 1,
                Annulus = 2,
                Column = 3,
                CementInner = 4,
                CementOuter = 5
            };
        };
        struct Thickness : public SomeProperty
        {
        };
        struct InnerRadius : public SomeProperty
        {
        };
        struct OuterRadius : public SomeProperty
        {
        };
        struct Depth : public SomeProperty
        {
        };

        struct Flow : public StationaryPhaseProperties
        {
            using StationaryPhaseProperties::StationaryPhaseProperties;
        };

        struct Tube : public StationaryPhaseProperties
        {
            using StationaryPhaseProperties::StationaryPhaseProperties;
        };

        struct Annulus : public StationaryPhaseProperties
        {
            using StationaryPhaseProperties::StationaryPhaseProperties;
        };

        struct Column : public StationaryPhaseProperties
        {
            using StationaryPhaseProperties::StationaryPhaseProperties;
        };

        struct Cement : public StationaryPhaseProperties
        {
            using StationaryPhaseProperties::StationaryPhaseProperties;
        };

        struct Ring
            : public StationaryPhaseProperties
        {
            template <typename PhaseProperties_t>
            Ring(const PhaseProperties_t &props,
                 Thickness thickness,
                 InnerRadius inner_radius,
                 Depth depth)
                : StationaryPhaseProperties{props},
                  thickness{thickness},
                  inner_radius{inner_radius},
                  outer_radius{inner_radius + thickness},
                  depth{depth}
            {
            }

            const RealType thickness, depth;
            const RealType inner_radius, outer_radius;
        };

        struct Isolation
        {
            Isolation(const std::vector<Ring> &completion)
                : sandwich{completion},
                  sandface_radius{completion.back().outer_radius},
                  flow_radius{completion.front().outer_radius}
            {
            }

            const std::vector<Ring> sandwich;
            const RealType sandface_radius, flow_radius;
        };

    } // Completion

} // GPN