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

            const auto &flow() const { return front(); }

            const auto &operator[](auto i) const
            {
                return sandwich[i];
            }

            const RealType casing_volumetric_heat_capacity() const
            {
                // exclude "flow" at "i = 0" from summation!
                RealType C{0.0};
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    C += m.linear_heat_capacity;
                }
                C /= casing_area();
                return C;
            }

            const RealType cement_volumetric_heat_capacity() const
            {
                // exclude "flow" at "i = 0" from summation!
                RealType C{0.0};
                for (ptrdiff_t i{MaterialType::Cement}; i <= MaterialType::Cement; ++i)
                {
                    const auto &m = sandwich[i];
                    C += m.linear_heat_capacity;
                }
                C /= cement_area();
                return C;
            }

            const RealType integral_casing_radial_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                RealType L{0.0};
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Column; ++i)
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

            const RealType radial_node_position() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                const auto &column{sandwich[MaterialType::Column]};

                const RealType T{T_avg()};

                const RealType R{tube.inner_radius * std::exp(T)};
                assert(R > tube.inner_radius);
                if (R < tube.outer_radius)
                    return R;
                else
                {
                    const RealType R{
                        annulus.inner_radius *
                        std::exp(
                            (T - std::log(tube.outer_radius / tube.inner_radius)) /
                            (tube.heat_conductivity / annulus.heat_conductivity))};
                    assert(R > annulus.inner_radius);
                    if (R < column.inner_radius)
                        return R;
                    else
                    {
                        const RealType R{
                            column.inner_radius *
                            std::exp(
                                (T - std::log(tube.outer_radius / tube.inner_radius) -
                                 (tube.heat_conductivity / annulus.heat_conductivity) *
                                     std::log(annulus.outer_radius / annulus.inner_radius)) /
                                (tube.heat_conductivity / column.heat_conductivity))};
                        assert(R > column.inner_radius);
                        assert(R < column.outer_radius);
                        return R;
                    }
                }
            }

            const RealType zeta_0(const RealType r1) const
            {
                const auto &Tube{sandwich[MaterialType::Tube]};
                const auto &Annulus{sandwich[MaterialType::Annulus]};
                const auto &Column{sandwich[MaterialType::Column]};

                return r1 < Tube.outer_radius
                           ? std::log(r1 / Tube.inner_radius) / Tube.heat_conductivity
                       : (r1 < Column.inner_radius)
                           ? 1 / Tube.radial_heat_conductivity + std::log(r1 / Tube.outer_radius) / Annulus.heat_conductivity
                           : 1 / Tube.radial_heat_conductivity + 1 / Annulus.radial_heat_conductivity + std::log(r1 / Column.inner_radius) / Column.heat_conductivity;
            }

            const RealType zeta_02(const RealType r2) const
            {
                const auto &Column{sandwich[MaterialType::Column]};
                const auto &Sandface{sandwich[MaterialType::Cement]};
                return 1 / integral_casing_radial_heat_conductivity() +
                                       std::log(r2 / Column.outer_radius) / Sandface.heat_conductivity;
            }

        private:
            const size_t size() const
            {
                return sandwich.size();
            }

            RealType I_tube() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                return tube.outer_radius * tube.outer_radius *
                           std::log(tube.outer_radius / tube.inner_radius) -
                       (tube.outer_radius - tube.inner_radius) * (tube.outer_radius + tube.inner_radius) / 2.0;
            }

            RealType I_annulus() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                return (annulus.area() / std::numbers::pi) *
                           std::log(tube.outer_radius / tube.inner_radius) +
                       tube.heat_conductivity / annulus.heat_conductivity *
                           (annulus.outer_radius * annulus.outer_radius *
                                std::log(annulus.outer_radius / tube.outer_radius) -
                            annulus.area() / 2.0 / std::numbers::pi);
            }

            RealType I_column() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                const auto &column{sandwich[MaterialType::Column]};
                return (
                           std::log(tube.outer_radius / tube.inner_radius) +
                           tube.heat_conductivity / annulus.heat_conductivity *
                               std::log(annulus.outer_radius / annulus.inner_radius)) *
                           (column.area() / std::numbers::pi) +
                       tube.heat_conductivity / column.heat_conductivity *
                           (column.outer_radius * column.outer_radius * std::log(column.outer_radius / column.inner_radius) -
                            column.area() / 2.0 / std::numbers::pi);
            }

            RealType T_avg() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                const auto &column{sandwich[MaterialType::Column]};

                return (tube.volumetric_heat_capacity * I_tube() +
                        annulus.volumetric_heat_capacity * I_annulus() +
                        column.volumetric_heat_capacity * I_column()) /
                       (casing_area() * casing_volumetric_heat_capacity());
            }
        };

    } // Completion

} // GPN

#pragma warning(pop)