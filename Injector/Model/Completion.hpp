#pragma once

#include <numbers>
#include <algorithm>

#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Properties/Logs.hpp>

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
                Cement = 4,
                Size = 5
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

        struct SomePropertyVar
        {
            operator Eigen::ArrayX<RealType>() const { return value.data; }
            Logs::StepProperty value;
        };

        auto operator*(const Logs::StepProperty &lhs, const Logs::StepProperty &rhs)
        {
            return Logs::StepProperty{StepPropertyContainer{lhs.data * rhs.data}};
        }
        auto operator-(const Logs::StepProperty &lhs, const Logs::StepProperty &rhs)
        {
            return Logs::StepProperty{StepPropertyContainer{lhs.data - rhs.data}};
        }
        auto operator+(const Logs::StepProperty &lhs, const Logs::StepProperty &rhs)
        {
            return Logs::StepProperty{StepPropertyContainer{lhs.data + rhs.data}};
        }
        auto operator*(RealType v, const Logs::StepProperty &rhs)
        {
            return Logs::StepProperty{StepPropertyContainer{v * rhs.data}};
        }
        auto operator*(const Logs::StepProperty &lhs, RealType v)
        {
            return v * lhs;
        }

        struct Density : public SomePropertyVar
        {
        };
        struct SpecificHeatCapacity : public SomePropertyVar
        {
        };
        struct VolumetricHeatCapacity : public SomePropertyVar
        {
        };
        struct HeatConductivity : public SomePropertyVar
        {
        };

        struct VarThickness : public SomePropertyVar
        {
        };
        struct VarInnerRadius : public SomePropertyVar
        {
        };
        struct VarDepth : public SomePropertyVar
        {
        };

        struct VariableStationaryPhaseProperties
        {
            VariableStationaryPhaseProperties(
                const VariableStationaryPhaseProperties &) = default;

            VariableStationaryPhaseProperties(
                const Density &density,
                const SpecificHeatCapacity &mass_heat_capacity,
                const HeatConductivity &heat_conductivity) noexcept
                : density{density},
                  mass_heat_capacity{mass_heat_capacity},
                  heat_conductivity{heat_conductivity},
                  volumetric_heat_capacity{
                      (mass_heat_capacity.value *
                      density.value).data}
            {
            }

            const Eigen::ArrayX<RealType> density;
            const Eigen::ArrayX<RealType> mass_heat_capacity;
            const Eigen::ArrayX<RealType> heat_conductivity;
            const Eigen::ArrayX<RealType> volumetric_heat_capacity;

            //    const Container_t<RealType> depth_grid;
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

        struct VarCement : public VariableStationaryPhaseProperties
        {
            using VariableStationaryPhaseProperties::VariableStationaryPhaseProperties;
        };

        struct VarColumn : public VariableStationaryPhaseProperties
        {
            using VariableStationaryPhaseProperties::VariableStationaryPhaseProperties;
        };
        
        struct VarAnnulus : public VariableStationaryPhaseProperties
        {
            using VariableStationaryPhaseProperties::VariableStationaryPhaseProperties;
        };
        
        struct VarTube : public VariableStationaryPhaseProperties
        {
            using VariableStationaryPhaseProperties::VariableStationaryPhaseProperties;
        };
        
        struct VarFlow : public VariableStationaryPhaseProperties
        {
            using VariableStationaryPhaseProperties::VariableStationaryPhaseProperties;
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

        struct VarRing
            : public VariableStationaryPhaseProperties
        {
            VarRing(const VariableStationaryPhaseProperties &props,
                    const VarThickness &thickness,
                    const VarInnerRadius &inner_radius,
                    const VarDepth &depth)
                : VariableStationaryPhaseProperties{props},
                  thickness{thickness},
                  inner_radius{inner_radius},
                  outer_radius{(inner_radius.value + thickness.value).data},
                  depth{depth},
                  linear_heat_capacity{
                        linear_heat_capacity_calc(
                            inner_radius, thickness, props)
                  },
                  radial_heat_conductivity{
                        radial_heat_conductivity_calc(
                            inner_radius, thickness, props)
                  },
                  integral_vertical_heat_conductivity{
                        props.heat_conductivity * 
                        std::numbers::pi *
                       ( thickness.value * (thickness.value + 2.0 * inner_radius.value)).data
                  }
            {
                   assert(assertion());
            }

            const Eigen::ArrayX<RealType>
                thickness, depth, inner_radius, outer_radius;

            // c
            const Eigen::ArrayX<RealType> linear_heat_capacity;
            // lambda/log(r_o/r_i)
            const Eigen::ArrayX<RealType> radial_heat_conductivity;
            const Eigen::ArrayX<RealType> integral_vertical_heat_conductivity;

            const Eigen::ArrayX<RealType> area() const
            {
                return std::numbers::pi *
                       (outer_radius - inner_radius) *
                       (outer_radius + inner_radius);
            }

        private:
            static Eigen::ArrayX<RealType>
            linear_heat_capacity_calc(
                VarInnerRadius inner_radius,
                VarThickness thickness,
                const VariableStationaryPhaseProperties &props)
            {
                return (std::numbers::pi *
                       thickness.value * (2.0 * inner_radius.value + thickness.value) *
                       props.volumetric_heat_capacity).data;
            }

            static Eigen::ArrayX<RealType>
            radial_heat_conductivity_calc(
                VarInnerRadius inner_radius,
                VarThickness thickness,
                const VariableStationaryPhaseProperties &props)
            {
                Eigen::ArrayX<RealType> out(thickness.value.size());
                for(auto idx{0ll}; idx < out.size(); ++idx)

                if (thickness.value(idx) == 0.0)
                {
                    out[idx] = std::numeric_limits<RealType>::infinity();
                }
                else
                {
                    const RealType temp{std::log(1.0 + (RealType)thickness.value(idx) / (RealType)inner_radius.value(idx))};
                    out[idx] = props.heat_conductivity(idx) / temp;
                }

                return out;
            }

            bool assertion()
            {
                const auto temp{
                    (outer_radius - inner_radius - thickness)};

                return std::all_of(
                    temp.cbegin(), temp.cend(),
                    [](RealType v)
                    { return std::abs(v) < 1e-12; });
            }
        };

        struct FlowRing : public Ring
        {
            using Ring::Ring;
            FlowRing(const Ring &r) : Ring{r} {}
        };
        struct TubeRing : public Ring
        {
            using Ring::Ring;
            TubeRing(const Ring &r) : Ring{r} {}
        };
        struct AnnulusRing : public Ring
        {
            using Ring::Ring;
            AnnulusRing(const Ring &r) : Ring{r} {}
        };
        struct ColumnRing : public Ring
        {
            using Ring::Ring;
            ColumnRing(const Ring &r) : Ring{r} {}
        };
        struct CementRing : public Ring
        {
            using Ring::Ring;
            CementRing(const Ring &r) : Ring{r} {}
        };

        struct Casing
        {
            Casing(const std::vector<Ring> &completion)
                : sandwich{completion},
                  sandface_radius{completion.back().outer_radius},
                  flow_radius{completion.front().outer_radius},
                  column_outer_radius{completion[MaterialType::Column].outer_radius}
            {
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Cement; ++i)
                {
                    assert(std::abs(completion[i].inner_radius - completion[i - 1ll].outer_radius) < 1e-12);
                }
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
            const RealType sandface_radius, column_outer_radius, flow_radius //, thickness
                ;

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

        struct ExtrudedRing : public StationaryPhaseProperties
        {
            template <typename Container_t, typename PhaseProperties_t>
            ExtrudedRing(
                const PhaseProperties_t &props,
                const Container_t &inner_radius,
                const Container_t &thickness)
                : StationaryPhaseProperties{props},
                  thickness{thickness},
                  inner_radius{inner_radius},
                  outer_radius{inner_radius + thickness},
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
                for (auto id{0ll}; id < thickness.size(); ++id)
                    assert(std::abs(outer_radius(id) - inner_radius(id) - thickness(id)) < 1e-12);
            }

            const Eigen::ArrayX<RealType>
                thickness, depth,
                inner_radius, outer_radius;

            // c
            const Eigen::ArrayX<RealType> linear_heat_capacity;
            // lambda/log(r_o/r_i)
            const Eigen::ArrayX<RealType> radial_heat_conductivity;
            const Eigen::ArrayX<RealType> integral_vertical_heat_conductivity;

            const Eigen::ArrayX<RealType> area() const
            {
                return std::numbers::pi *
                       (outer_radius - inner_radius) *
                       (outer_radius + inner_radius);
            }

        private:
            template <typename PhaseProperties_t>
            static auto
            linear_heat_capacity_calc(
                const auto &inner_radius,
                const auto &thickness,
                const PhaseProperties_t &props)
            {
                return std::numbers::pi *
                       thickness * (2 * inner_radius + thickness) *
                       props.volumetric_heat_capacity;
            }

            template <typename PhaseProperties_t>
            static auto
            radial_heat_conductivity_calc(
                const auto &inner_radius,
                const auto &thickness,
                const PhaseProperties_t &props)
            {
                const auto temp{(1.0 + thickness / inner_radius).log()};
                Eigen::ArrayX<RealType> out{props.heat_conductivity / temp};
                for (auto id{0ll}; id < out.size(); ++id)
                {
                    if (thickness(id) == 0.0)
                        out(id) = std::numeric_limits<RealType>::infinity();
                }

                return out;
            }
        };

        struct ExtrudedCasing
        {
            ExtrudedCasing(const std::vector<ExtrudedRing> &completion)
                : sandwich{completion},
                  sandface_radius{completion.back().outer_radius},
                  flow_radius{completion.front().outer_radius},
                  column_outer_radius{completion[MaterialType::Column].outer_radius}
            {
            }

            const std::vector<ExtrudedRing> sandwich;
            const Eigen::ArrayX<RealType> sandface_radius, column_outer_radius;
            const Eigen::ArrayX<RealType> flow_radius;

            const auto &back() const { return sandwich.back(); }
            const auto &front() const { return sandwich.front(); }

            const auto &flow() const { return front(); }

            const auto &operator[](auto i) const
            {
                return sandwich[i];
            }

            const Eigen::ArrayX<RealType> cement_area() const
            {
                const auto out{std::numbers::pi *
                               (sandwich[MaterialType::Cement].thickness) *
                               (sandwich[MaterialType::Cement].outer_radius + sandwich[MaterialType::Cement].inner_radius)};
                assert(std::all_of(out.cbegin(), out.cend(), [](const auto v)
                                   { return v > 0.0; }));
                return out;
            }

            const Eigen::ArrayX<RealType> casing_area() const
            {
                const auto out{std::numbers::pi *
                               (sandwich[MaterialType::Column].outer_radius - sandwich[MaterialType::Tube].inner_radius) *
                               (sandwich[MaterialType::Column].outer_radius + sandwich[MaterialType::Tube].inner_radius)};
                assert(std::all_of(out.cbegin(), out.cend(), [](const auto v)
                                   { return v > 0.0; }));
                return out;
            }

            const Eigen::ArrayX<RealType> casing_volumetric_heat_capacity() const
            {
                // exclude "flow" at "i = 0" from summation!
                auto C{sandwich[MaterialType::Tube].linear_heat_capacity};
                for (ptrdiff_t i{MaterialType::Tube + 1ll}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    C += m.linear_heat_capacity;
                }
                C /= casing_area();
                return C;
            }

            const Eigen::ArrayX<RealType> cement_volumetric_heat_capacity() const
            {
                return sandwich[MaterialType::Cement].linear_heat_capacity / cement_area();
            }

            const Eigen::ArrayX<RealType> integral_casing_radial_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                Eigen::ArrayX<RealType> L{1.0 / sandwich[MaterialType::Tube].radial_heat_conductivity};
                for (ptrdiff_t i{MaterialType::Tube + 1ll}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    L += 1.0 / m.radial_heat_conductivity;
                }
                assert(std::all_of(L.cbegin(), L.cend(), [](const RealType v)
                                   { return (v > 0.0) && !std::isnan(v); }));
                return 1.0 / L;
            }

            const Eigen::ArrayX<RealType> integral_vertical_casing_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                Eigen::ArrayX<RealType> L{sandwich[MaterialType::Tube].integral_vertical_heat_conductivity};
                for (ptrdiff_t i{MaterialType::Tube + 1ll}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    L += m.integral_vertical_heat_conductivity;
                }
                return L / casing_area();
            }

            const Eigen::ArrayX<RealType> integral_vertical_cement_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                Eigen::ArrayX<RealType> L{sandwich[MaterialType::Cement].integral_vertical_heat_conductivity};
                return L / cement_area();
            }

            const Eigen::ArrayX<RealType> area() const
            {
                return cement_area() + casing_area();
            }

            const Eigen::ArrayX<RealType> I_tube() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                return tube.outer_radius * tube.outer_radius *
                           log(tube.outer_radius / tube.inner_radius) -
                       (tube.outer_radius - tube.inner_radius) * (tube.outer_radius + tube.inner_radius) / 2.0;
            }

            const Eigen::ArrayX<RealType> I_annulus() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                return (annulus.area() / std::numbers::pi) *
                           log(tube.outer_radius / tube.inner_radius) +
                       tube.heat_conductivity / annulus.heat_conductivity *
                           (annulus.outer_radius * annulus.outer_radius *
                                log(annulus.outer_radius / tube.outer_radius) -
                            annulus.area() / 2.0 / std::numbers::pi);
            }

            const Eigen::ArrayX<RealType> I_column() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                const auto &column{sandwich[MaterialType::Column]};
                return (
                           log(tube.outer_radius / tube.inner_radius) +
                           tube.heat_conductivity / annulus.heat_conductivity *
                               log(annulus.outer_radius / annulus.inner_radius)) *
                           (column.area() / std::numbers::pi) +
                       tube.heat_conductivity / column.heat_conductivity *
                           (column.outer_radius * column.outer_radius * log(column.outer_radius / column.inner_radius) -
                            column.area() / 2.0 / std::numbers::pi);
            }

            const Eigen::ArrayX<RealType> T_avg() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                const auto &column{sandwich[MaterialType::Column]};

                return (tube.volumetric_heat_capacity * I_tube() +
                        annulus.volumetric_heat_capacity * I_annulus() +
                        column.volumetric_heat_capacity * I_column()) /
                       (casing_area() * casing_volumetric_heat_capacity());
            }

            const Eigen::ArrayX<RealType> radial_node_position() const
            {
                const auto &tube{sandwich[MaterialType::Tube]};
                const auto &annulus{sandwich[MaterialType::Annulus]};
                const auto &column{sandwich[MaterialType::Column]};

                const auto T{T_avg()};

                Eigen::ArrayX<RealType> out(T.size());

                for (auto id{0ll}; id < T.size(); ++id)
                {
                    const auto R{tube.inner_radius(id) * std::exp(T(id))};
                    assert(R > tube.inner_radius(id));
                    if (R < tube.outer_radius(id))
                        out(id) = R;
                    else
                    {
                        const RealType R{
                            annulus.inner_radius(id) *
                            exp(
                                (T(id) - std::log(tube.outer_radius(id) / tube.inner_radius(id))) /
                                (tube.heat_conductivity / annulus.heat_conductivity))};
                        assert(R > annulus.inner_radius(id));
                        if (R < column.inner_radius(id))
                            out(id) = R;
                        else
                        {
                            const RealType R{
                                column.inner_radius(id) *
                                std::exp(
                                    (T(id) - std::log(tube.outer_radius(id) / tube.inner_radius(id)) -
                                     (tube.heat_conductivity / annulus.heat_conductivity) *
                                         std::log(annulus.outer_radius(id) / annulus.inner_radius(id))) /
                                    (tube.heat_conductivity / column.heat_conductivity))};
                            assert(R > column.inner_radius(id));
                            assert(R < column.outer_radius(id));
                            out(id) = R;
                        }
                    }
                }
                return out;
            }

            const Eigen::ArrayX<RealType> zeta_0(const auto r1) const
            {
                const auto &Tube{sandwich[MaterialType::Tube]};
                const auto &Annulus{sandwich[MaterialType::Annulus]};
                const auto &Column{sandwich[MaterialType::Column]};

                Eigen::ArrayX<RealType> out(r1.size());

                for (auto id{0ll}; id < r1.size(); ++id)
                    out(id) = r1(id) < Tube.outer_radius(id)
                                  ? log(r1(id) / Tube.inner_radius(id)) / Tube.heat_conductivity
                              : (r1(id) < Column.inner_radius(id))
                                  ? 1 / Tube.radial_heat_conductivity(id) + log(r1(id) / Tube.outer_radius(id)) / Annulus.heat_conductivity
                                  : 1 / Tube.radial_heat_conductivity(id) + 1 / Annulus.radial_heat_conductivity(id) + log(r1(id) / Column.inner_radius(id)) / Column.heat_conductivity;
                return out;
            }

            const Eigen::ArrayX<RealType> zeta_02(const RealType r2) const
            {
                const auto &Column{sandwich[MaterialType::Column]};
                const auto &Sandface{sandwich[MaterialType::Cement]};
                return 1 / integral_casing_radial_heat_conductivity() +
                       log(r2 / Column.outer_radius) / Sandface.heat_conductivity;
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