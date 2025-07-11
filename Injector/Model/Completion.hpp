#pragma once

#include <numbers>

#include <Injector/Model/Phases/PhaseProperties.hpp>
#pragma warning(push)
#pragma warning(disable : 4723)
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
                Cement = 4
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
                  depth{depth},
                  linear_heat_capacity{
                      linear_heat_capacity_calc(
                          inner_radius, thickness, props)},
                  radial_heat_conductivity{
                      radial_heat_conductivity_calc(
                          inner_radius, thickness, props)},
                  integral_vertical_heat_conductivity{
                      props.heat_conductivity * std::numbers::pi *
                      thickness * (thickness + 2.0 * inner_radius)}
            {
                assert(std::abs(outer_radius - inner_radius - thickness) < 1e-12);
            }

            const RealType thickness, depth;
            const RealType inner_radius, outer_radius;

            // c
            const RealType linear_heat_capacity;
            // lambda/log(r_o/r_i)
            const RealType radial_heat_conductivity;
            const RealType integral_vertical_heat_conductivity;

            const RealType area() const
            {
                return std::numbers::pi *
                       (outer_radius - inner_radius) *
                       (outer_radius + inner_radius);
            }

        private:
            template <typename PhaseProperties_t>
            static RealType
            linear_heat_capacity_calc(
                InnerRadius inner_radius,
                Thickness thickness,
                const PhaseProperties_t &props)
            {
                return std::numbers::pi *
                       thickness * (2 * inner_radius + thickness) *
                       props.volumetric_heat_capacity;
            }

            template <typename PhaseProperties_t>
            static RealType
            radial_heat_conductivity_calc(
                InnerRadius inner_radius,
                Thickness thickness,
                const PhaseProperties_t &props)
            {
                if (thickness == 0.0)
                {
                    return std::numeric_limits<RealType>::infinity();
                }
                else
                {
                    const RealType temp{std::log((1.0 + (RealType)thickness / (RealType)inner_radius))};
                    return props.heat_conductivity / temp;
                }
            }
        };

        struct Casing
        {
            Casing(const std::vector<Ring> &completion)
                : sandwich{completion},
                  sandface_radius{completion.back().outer_radius},
                  flow_radius{completion.front().outer_radius},
                  column_outer_radius{completion[MaterialType::Column].outer_radius},
                  thickness{completion.back().outer_radius - completion.front().outer_radius}
            {
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Cement; ++i)
                {
                    assert(std::abs(completion[i].inner_radius - completion[i - 1ll].outer_radius) < 1e-12);
                }

                assert(std::abs(thickness - sandface_radius + flow_radius) < 1e-12);
            }

            const auto &back() const { return sandwich.back(); }
            const auto &front() const { return sandwich.front(); }

            const auto &operator[](auto i) const
            {
                return sandwich[i];
            }

            const RealType volumetric_heat_capacity() const
            {
                // exclude "flow" at "i = 0" from summation!
                RealType C{0.0};
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Cement; ++i)
                {
                    const auto &m = sandwich[i];
                    C += m.linear_heat_capacity;
                }
                C /= area();
                return C;
            }

            const RealType integral_inner_radial_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                RealType L{0.0};
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Cement; ++i)
                {
                    const auto &m = sandwich[i];
                    L += 1.0 / m.radial_heat_conductivity;
                }
                assert(L > 0.0);
                assert(!std::isnan(L));
                return 1.0 / L;
            }

            const RealType integral_vertical_casing_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                RealType L{0.0};
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    L += m.integral_vertical_heat_conductivity;
                }
                return L / casing_area();
            }

            const RealType integral_vertical_cement_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                RealType L{0.0};
                for (ptrdiff_t i{MaterialType::Cement}; i <= MaterialType::Cement; ++i)
                {
                    const auto &m = sandwich[i];
                    L += m.integral_vertical_heat_conductivity;
                }
                return L / cement_area();
            }

            const std::vector<Ring> sandwich;
            const RealType sandface_radius, column_outer_radius, flow_radius, thickness;

            const RealType area() const
            {
                const auto out{cement_area() + casing_area()};
                assert(out > 0.0);
                return out;
            }

            const RealType cement_area() const
            {
                const auto out{std::numbers::pi *
                               (sandwich[MaterialType::Cement].thickness) *
                               (sandwich[MaterialType::Cement].outer_radius + sandwich[MaterialType::Cement].inner_radius)};
                assert(out > 0.0);
                return out;
            }
            
            const RealType casing_area() const
            {
                const auto out{std::numbers::pi *
                               (sandwich[MaterialType::Column].outer_radius - sandwich[MaterialType::Tube].inner_radius) *
                               (sandwich[MaterialType::Column].outer_radius + sandwich[MaterialType::Tube].inner_radius)};
                assert(out > 0.0);
                return out;
            }

        private:
            const size_t size() const
            {
                return sandwich.size();
            }
        };

    } // Completion

} // GPN

#pragma warning(pop)