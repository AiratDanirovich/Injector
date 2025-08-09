#pragma once

#include <numbers>
#include <algorithm>
#include <vector>
#include <exception>

#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Properties/Logs.hpp>

#include <tests/includes/print_vector.hpp>
#include <tests/includes/accumulate_steps.hpp>

#pragma warning(push)
#pragma warning(disable : 4723)
namespace GPN
{
    namespace Completion
    {
        struct MaterialType
        {
            enum material_type
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

        /// @brief Intervals of depth of const-value properties
        struct VarDepth : public SomePropertyVar
        {
        };

        struct VariableStationaryPhaseProperties
        {
            VariableStationaryPhaseProperties(
                const VariableStationaryPhaseProperties &) = default;

            VariableStationaryPhaseProperties(
                const Density &density,
                const SpecificHeatCapacity &specific_heat_capacity,
                const HeatConductivity &heat_conductivity) noexcept
                : density{density},
                  specific_heat_capacity{specific_heat_capacity},
                  heat_conductivity{heat_conductivity},
                  volumetric_heat_capacity{
                      (specific_heat_capacity.value *
                       density.value)
                          .data}
            {
            }

            /// @brief Generate const-value phase properties from a single set of parameters
            /// @param density
            /// @param specific_heat_capacity
            /// @param heat_conductivity
            VariableStationaryPhaseProperties(
                const GPN::Density &density,
                const GPN::SpecificHeatCapacity &specific_heat_capacity,
                const GPN::HeatConductivity &heat_conductivity)
                : VariableStationaryPhaseProperties{
                      Density{std::vector<RealType>{density}},
                      SpecificHeatCapacity{std::vector<RealType>{specific_heat_capacity}},
                      HeatConductivity{std::vector<RealType>{heat_conductivity}}}
            {
            }

            const Eigen::ArrayX<RealType> density;
            const Eigen::ArrayX<RealType> specific_heat_capacity;
            const Eigen::ArrayX<RealType> heat_conductivity;
            const Eigen::ArrayX<RealType> volumetric_heat_capacity;
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
            using value_type = RealType;

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
                  max_depth{depth},
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

            const RealType thickness, depth, max_depth;
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

        struct VarRingSimple
        {
            using value_type = Eigen::ArrayX<RealType>;

            VarRingSimple(const VariableStationaryPhaseProperties &props,
                          const VarThickness &thickness,
                          const VarDepth &depth_intervals)
                : props{props},
                  thickness{thickness},
                  depth_stencils{depth_stencils_calc(depth_intervals)},
                  max_depth{depth_intervals.value.data.sum()}
            {
            }

            const auto &density() const
            {
                return props.density;
            }
            const auto &specific_heat_capacity() const
            {
                return props.specific_heat_capacity;
            }
            const auto &volumetric_heat_capacity() const
            {
                return props.volumetric_heat_capacity;
            }
            const auto &heat_conductivity() const
            {
                return props.heat_conductivity;
            }

            const VariableStationaryPhaseProperties props;

            const value_type
                thickness,
                depth_stencils;
            const RealType max_depth;

        private:
            static value_type
            depth_stencils_calc(const VarDepth &depth)
            {
                const auto &data{depth.value.data};
                value_type out(1ll + data.size());
                out(0ll) = 0.0;
                for (auto i{1ll}; i < out.size(); ++i)
                    out(i) = out(i - 1ll) + data(i - 1ll);
                return out;
            }
        };

        struct VarRing : public VarRingSimple
        {
            using value_type = Eigen::ArrayX<RealType>;

            VarRing(const VarRingSimple &props,
                    const VarInnerRadius &inner_radius)
                : VarRingSimple{props},
                  inner_radius{inner_radius},
                  outer_radius{(value_type)inner_radius + (value_type)thickness},
                  linear_heat_capacity{
                      linear_heat_capacity_calc(
                          this->inner_radius, this->thickness, props.volumetric_heat_capacity())},
                  radial_heat_conductivity{
                      radial_heat_conductivity_calc(
                          this->inner_radius, this->thickness, props.heat_conductivity())},
                  integral_vertical_heat_conductivity{
                      props.heat_conductivity() *
                      std::numbers::pi *
                      (value_type)thickness * ((value_type)thickness + 2.0 * (value_type)inner_radius)}
            {
                assert(assertion());
            }

            VarRing(const VariableStationaryPhaseProperties &props,
                    const VarThickness &thickness,
                    const VarInnerRadius &inner_radius,
                    const VarDepth &depth_intervals)
                : VarRing{
                      VarRingSimple{props, thickness, depth_intervals},
                      inner_radius}
            {
            }

            // static value_type populate_array(const value_type &property, const auto &depth_stencils, const auto &grid_z)
            // {
            //     value_type out(grid_z.mesh_size());
            //     auto id{0ll};
            //     for (auto depth_id{0ll}; id < grid_z.mesh_size(); ++id)
            //     {
            //         if (depth_stencils(depth_id) < grid_z.mesh_nodes(id))
            //             ++depth_id;
            //         if (depth_stencils(depth_id) < grid_z.mesh_nodes(id))
            //             throw std::logic_error(
            //                 "Z-grid is too coarce! "
            //                 "The casing properties vary multile times within a single z-step.");
            //         out(id) = property(depth_id);
            //     }
            //     return out;
            // }

            const value_type
                inner_radius,
                outer_radius;

            // c
            const value_type linear_heat_capacity;
            // lambda/log(r_o/r_i)
            const value_type radial_heat_conductivity;
            const value_type integral_vertical_heat_conductivity;

            const value_type area() const
            {
                return std::numbers::pi *
                       (outer_radius - inner_radius) *
                       (outer_radius + inner_radius);
            }

        private:
            static value_type
            depth_stencils_calc(const VarDepth &depth)
            {
                const auto &data{depth.value.data};
                value_type out(1ll + data.size());
                out(0ll) = 0.0;
                for (auto i{1ll}; i < out.size(); ++i)
                    out(i) = out(i - 1ll) + data(i - 1ll);
                return out;
            }

            static value_type
            linear_heat_capacity_calc(
                const value_type &inner_radius,
                const value_type &thickness,
                const value_type &volumetric_heat_capacity)
            {
                return (std::numbers::pi *
                        thickness * (2.0 * inner_radius + thickness) *
                        volumetric_heat_capacity);
            }

            static value_type
            radial_heat_conductivity_calc(
                const value_type &inner_radius,
                const value_type &thickness,
                const value_type &heat_conductivity)
            {
                Eigen::ArrayX<RealType> out(thickness.size());
                for (auto idx{0ll}; idx < out.size(); ++idx)

                    if (thickness(idx) == 0.0)
                    {
                        out[idx] = std::numeric_limits<RealType>::infinity();
                    }
                    else
                    {
                        const RealType temp{std::log(1.0 + (RealType)thickness(idx) / (RealType)inner_radius(idx))};
                        out[idx] = heat_conductivity(idx) / temp;
                    }

                return out;
            }

            bool assertion()
            {
                const value_type temp{
                    (outer_radius - inner_radius - thickness)};

                return std::all_of(
                    temp.cbegin(), temp.cend(),
                    [](RealType v)
                    { return std::abs(v) < 1e-12; });
            }
        };

        struct FactoryVarRing
        {
            static VarRing create_column_ring(
                const VarRing &column,
                const auto &grid_z)
            {
                using namespace std;

                vector<RealType>
                    density, specific_heat_capacity, heat_conductivity, outer_radius, inner_radius;
                density.reserve(grid_z.mesh_size());
                specific_heat_capacity.reserve(grid_z.mesh_size());
                heat_conductivity.reserve(grid_z.mesh_size());
                outer_radius.reserve(grid_z.mesh_size());
                inner_radius.reserve(grid_z.mesh_size());

                for (auto id{0ll}, depth_id{1ll}; id < grid_z.mesh_size(); ++id)
                {
                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        ++depth_id;
                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        throw std::logic_error(
                            "Column physical properties discretization is "
                            "higher than the discretization of z-axis grid.");

                    density.push_back(
                        column.density()(depth_id - 1ll));
                    specific_heat_capacity.push_back(
                        column.specific_heat_capacity()(depth_id - 1ll));
                    heat_conductivity.push_back(
                        column.heat_conductivity()(depth_id - 1ll));

                    outer_radius.push_back(
                        column.outer_radius(depth_id - 1ll));
                    inner_radius.push_back(
                        column.inner_radius(depth_id - 1ll));
                }

                for (auto id{0ull}; id < outer_radius.size(); ++id)
                    assert(outer_radius[id] > inner_radius[id]);

                vector<RealType> thickness(outer_radius.size());
                std::transform(
                    outer_radius.cbegin(), outer_radius.cend(),
                    inner_radius.cbegin(), thickness.begin(),
                    [](auto a, auto b)
                    { return a - b; });

                return VarRing{
                    VarRingSimple{
                        VarColumn{
                            Completion::Density{density},
                            Completion::SpecificHeatCapacity{specific_heat_capacity},
                            Completion::HeatConductivity{heat_conductivity}},
                        Completion::VarThickness{thickness},
                        Completion::VarDepth{grid_z.dual_steps}},
                    Completion::VarInnerRadius{inner_radius}};
            }

            static VarRing create_tube_ring(
                const VarRing &tube,
                const VarRing &column,
                const auto &grid_z)
            {
                using namespace std;
                
                assert(grid_z.mesh_back() > tube.max_depth);

                vector<RealType>
                    density, specific_heat_capacity, heat_conductivity,
                    outer_radius, inner_radius;
                density.reserve(grid_z.mesh_size());
                specific_heat_capacity.reserve(grid_z.mesh_size());
                heat_conductivity.reserve(grid_z.mesh_size());
                outer_radius.reserve(grid_z.mesh_size());
                inner_radius.reserve(grid_z.mesh_size());

                auto id{0ll};
                for (
                    auto depth_id{1ll};
                    (id < grid_z.mesh_size()) &&
                    (grid_z.mesh_nodes(id) < tube.max_depth);
                    ++id)
                {
                    if (grid_z.mesh_nodes(id) > tube.depth_stencils(depth_id))
                        ++depth_id;
                    if (grid_z.mesh_nodes(id) > tube.depth_stencils(depth_id))
                        throw std::logic_error(
                            "Tube physical properties discretization is "
                            "higher than the discretization of z-axis grid.");

                    density.push_back(
                        tube.density()(depth_id - 1ll));
                    specific_heat_capacity.push_back(
                        tube.specific_heat_capacity()(depth_id - 1ll));
                    heat_conductivity.push_back(
                        tube.heat_conductivity()(depth_id - 1ll));

                    outer_radius.push_back(
                        tube.outer_radius(depth_id - 1ll));
                    inner_radius.push_back(
                        tube.inner_radius(depth_id - 1ll));
                }
                // extrapolate tube physical properties
                // below tube length
                // with the last set of (density; capacity; conductivity) annulus values
                auto depth_id{1ll};
                while (column.depth_stencils(depth_id) < grid_z.mesh_nodes(id))
                    ++depth_id;

                for (; id < grid_z.mesh_size(); ++id)
                {
                    density.push_back(
                        tube.density().tail(1ll)(0ll));
                    specific_heat_capacity.push_back(
                        tube.specific_heat_capacity().tail(1ll)(0ll));
                    heat_conductivity.push_back(
                        tube.heat_conductivity().tail(1ll)(0ll));
                        
                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        ++depth_id;
                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        throw std::logic_error(
                            "Column physical properties discretization is "
                            "higher than the discretization of z-axis grid.");

                    outer_radius.push_back(
                        column.inner_radius(depth_id - 1ll));
                    inner_radius.push_back(
                        column.inner_radius(depth_id - 1ll));
                }
                for (auto id{0ull}; id < outer_radius.size(); ++id)
                {
                    if (grid_z.mesh_nodes(id) < tube.max_depth)
                        assert(outer_radius[id] > inner_radius[id]);
                    else
                        assert(outer_radius[id] == inner_radius[id]);
                }

                vector<RealType> thickness(outer_radius.size());
                std::transform(
                    outer_radius.cbegin(), outer_radius.cend(),
                    inner_radius.cbegin(), thickness.begin(),
                    [](auto a, auto b)
                    { return a - b; });

                return VarRing{
                    VarRingSimple{
                        VarTube{
                            Completion::Density{density},
                            Completion::SpecificHeatCapacity{specific_heat_capacity},
                            Completion::HeatConductivity{heat_conductivity}},
                        Completion::VarThickness{thickness},
                        Completion::VarDepth{grid_z.dual_steps}},
                    Completion::VarInnerRadius{inner_radius}};
            }

            static VarRing create_flow_ring(
                const Flow &flow_ring_props,
                const VarRing &tube,
                const VarRing &column,
                const auto &grid_z)
            {
                using namespace std;

                vector<RealType>
                    density, specific_heat_capacity, heat_conductivity,
                    inner_radius, outer_radius;
                density.reserve(grid_z.mesh_size());
                specific_heat_capacity.reserve(grid_z.mesh_size());
                heat_conductivity.reserve(grid_z.mesh_size());
                inner_radius.reserve(grid_z.mesh_size());
                outer_radius.reserve(grid_z.mesh_size());

                assert(tube.depth_stencils.size() > 1ll);
                assert(column.depth_stencils.size() > 1ll);

                assert(grid_z.mesh_back() < column.max_depth);

                for (auto id{0ll}; id < grid_z.mesh_size(); ++id)
                {
                    density.push_back(flow_ring_props.density);
                    specific_heat_capacity.push_back(flow_ring_props.specific_heat_capacity);
                    heat_conductivity.push_back(flow_ring_props.heat_conductivity);
                    inner_radius.push_back(0.0);
                        outer_radius.push_back(tube.inner_radius(id));
                }

                for (auto id{0ull}; id < outer_radius.size(); ++id)
                    assert(outer_radius[id] > inner_radius[id]);

                vector<RealType> thickness(outer_radius.size());
                std::transform(
                    outer_radius.cbegin(), outer_radius.cend(),
                    inner_radius.cbegin(), thickness.begin(),
                    [](auto a, auto b)
                    { return a - b; });

                return VarRing{
                    VarRingSimple{
                        VarFlow{
                            Completion::Density{density},
                            Completion::SpecificHeatCapacity{specific_heat_capacity},
                            Completion::HeatConductivity{heat_conductivity}},
                        Completion::VarThickness{thickness},
                        Completion::VarDepth{grid_z.dual_steps}},
                    Completion::VarInnerRadius{inner_radius}};
            }

            static VarRing create_annulus_ring(
                const VarAnnulus &annulus_ring_props,
                const Eigen::ArrayX<RealType> &depth_intervals,
                const VarRing &tube,
                const VarRing &column,
                const auto &grid_z)
            {
                using namespace std;

                assert(tube.depth_stencils.size() > 1ll);
                assert(column.depth_stencils.size() > 1ll);

                const auto depth_stencils = accumulate_steps(depth_intervals);

                assert(depth_stencils.back() == tube.max_depth);

#pragma region SET-PHYSICAL-PROPERTIES
                vector<RealType>
                    density, specific_heat_capacity, heat_conductivity;
                density.reserve(grid_z.mesh_size());
                specific_heat_capacity.reserve(grid_z.mesh_size());
                heat_conductivity.reserve(grid_z.mesh_size());

                {
                    auto id{0ll};
                    for (
                        auto depth_id{1ll};
                        (id < grid_z.mesh_size()) &&
                        (grid_z.mesh_nodes(id) < depth_stencils.back());
                        ++id)
                    {
                        if (grid_z.mesh_nodes(id) > depth_stencils[depth_id])
                            ++depth_id;

                        if (grid_z.mesh_nodes(id) > depth_stencils[depth_id])
                            throw std::logic_error(
                                "Annulus physical properties discretization is "
                                "higher than the discretization of z-axis grid.");

                        density.push_back(
                            annulus_ring_props.density(depth_id - 1ll));
                        specific_heat_capacity.push_back(
                            annulus_ring_props.specific_heat_capacity(depth_id - 1ll));
                        heat_conductivity.push_back(
                            annulus_ring_props.heat_conductivity(depth_id - 1ll));
                    }

                    // extrapolate annulus physical properties
                    // below tube length
                    // with the last set of (density; capacity; conductivity) annulus values
                    for (; id < grid_z.mesh_size(); ++id)
                    {
                        density.push_back(
                            annulus_ring_props.density.tail(1ll)(0ll));
                        specific_heat_capacity.push_back(
                            annulus_ring_props.specific_heat_capacity.tail(1ll)(0ll));
                        heat_conductivity.push_back(
                            annulus_ring_props.heat_conductivity.tail(1ll)(0ll));
                    }
                    assert(density.size() == grid_z.mesh_size());
                    assert(specific_heat_capacity.size() == grid_z.mesh_size());
                    assert(heat_conductivity.size() == grid_z.mesh_size());
                }

                cout << "density:\n";
                print_vector(density);

                cout << "grid_z:\n";
                cout << grid_z.mesh_nodes.transpose() << endl;

                cout << "tube.depth_stencils:\n";
                cout << tube.depth_stencils.transpose() << endl;
#pragma endregion
#pragma region SET-INNER-RADIUS
                vector<RealType> inner_radius;
                inner_radius.reserve(grid_z.mesh_size());
                {
                    assert(grid_z.mesh_nodes.tail(1ll)(0ll) > tube.max_depth);
                    auto id{0ll};
                    for (
                        auto depth_id{1ll};
                        (id < grid_z.mesh_size()) &&
                        (grid_z.mesh_nodes(id) < depth_stencils.back());
                        ++id)
                    {
                        if (grid_z.mesh_nodes(id) > tube.depth_stencils(depth_id))
                            ++depth_id;

                        if (
                            (depth_id < tube.depth_stencils.size() - 1ll) &&
                            (grid_z.mesh_nodes(id) > tube.depth_stencils(depth_id)))
                            throw std::logic_error(
                                "Tube properties discretization is "
                                "higher than the discretization of z-axis grid.");

                        inner_radius.push_back(
                            tube.outer_radius(depth_id - 1ll));
                    }
                    // extrapolate annulus physical properties
                    // below tube length
                    // with the last set of (density; capacity; conductivity) annulus values
                    auto depth_id{1ll};
                    while (column.depth_stencils(depth_id) < grid_z.mesh_nodes(id))
                        ++depth_id;

                    inner_radius.push_back(
                        column.inner_radius(depth_id - 1ll));
                    ++id;

                    for (; id < grid_z.mesh_size(); ++id)
                    {
                        if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                            ++depth_id;

                        inner_radius.push_back(
                            column.inner_radius(depth_id - 1ll));
                    }
                    assert(inner_radius.size() == grid_z.mesh_size());
                }
#pragma endregion
#pragma region SET-OUTER-RADIUS
                vector<RealType> outer_radius;
                outer_radius.reserve(grid_z.mesh_size());
                assert(grid_z.mesh_nodes.tail(1ll)(0ll) < column.max_depth);
                for (auto id{0ll}, depth_id{1ll}; id < grid_z.mesh_size(); ++id)
                {
                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        ++depth_id;

                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        throw std::logic_error(
                            "Column properties discretization is "
                            "higher than the discretization of z-axis grid.");

                    outer_radius.push_back(
                        column.inner_radius(depth_id - 1ll));
                }
                assert(outer_radius.size() == grid_z.mesh_size());
#pragma endregion

                for (auto id{0ull}; id < outer_radius.size(); ++id)
                {
                    if (grid_z.mesh_nodes(id) < tube.max_depth)
                        assert(outer_radius[id] > inner_radius[id]);
                    else
                        assert(outer_radius[id] == inner_radius[id]);
                }

                vector<RealType> thickness(outer_radius.size());
                std::transform(
                    outer_radius.cbegin(), outer_radius.cend(),
                    inner_radius.cbegin(), thickness.begin(),
                    [](auto a, auto b)
                    { return a - b; });

                return VarRing{
                    VarRingSimple{
                        VarAnnulus{
                            Completion::Density{density},
                            Completion::SpecificHeatCapacity{specific_heat_capacity},
                            Completion::HeatConductivity{heat_conductivity}},
                        Completion::VarThickness{thickness},
                        Completion::VarDepth{grid_z.dual_steps}},
                    Completion::VarInnerRadius{inner_radius}};
            }

            static VarRing create_cement_ring(
                const VarCement &cement_ring_props,
                const Eigen::ArrayX<RealType> &depth_intervals,
                const VarRing &column,
                const auto &grid_z)
            {
                using namespace std;

                vector<RealType> depth_stencils(depth_intervals.size() + 1ull, 0.0);
                std::partial_sum(
                    depth_intervals.cbegin(),
                    depth_intervals.cend(),
                    depth_stencils.begin() + 1ull);
                assert(depth_stencils[0ull] == 0.0);
                for (auto i{1ull}; i < depth_stencils.size(); ++i)
                    assert(depth_stencils[i] == depth_stencils[i - 1] + depth_intervals(i - 1ll));

                assert(depth_stencils.back() == column.max_depth);

#pragma region SET-PHYSICAL-PROPERTIES
                vector<RealType>
                    density, specific_heat_capacity, heat_conductivity;
                density.reserve(grid_z.mesh_size());
                specific_heat_capacity.reserve(grid_z.mesh_size());
                heat_conductivity.reserve(grid_z.mesh_size());
                {
                }
#pragma endregion
#pragma region SET-OUTER-RADIUS
                vector<RealType> outer_radius;
                outer_radius.reserve(grid_z.mesh_size());
                assert(grid_z.mesh_nodes.tail(1ll)(0ll) < column.max_depth);
                for (auto id{0ll}, depth_id{1ll}; id < grid_z.mesh_size(); ++id)
                {
                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        ++depth_id;

                    if (grid_z.mesh_nodes(id) > column.depth_stencils(depth_id))
                        throw std::logic_error(
                            "Column properties discretization is "
                            "higher than the discretization of z-axis grid.");

                    outer_radius.push_back(
                        column.inner_radius(depth_id - 1ll));
                }
                assert(outer_radius.size() == grid_z.mesh_size());
#pragma endregion

                const auto &inner_radius{column.outer_radius};

                vector<RealType> thickness(outer_radius.size());
                std::transform(
                    outer_radius.cbegin(), outer_radius.cend(),
                    inner_radius.cbegin(), thickness.begin(),
                    [](auto a, auto b)
                    { return a - b; });

                return VarRing{
                    VarRingSimple{
                        VarCement{
                            Completion::Density{density},
                            Completion::SpecificHeatCapacity{specific_heat_capacity},
                            Completion::HeatConductivity{heat_conductivity}},
                        Completion::VarThickness{thickness},
                        Completion::VarDepth{grid_z.dual_steps}},
                    Completion::VarInnerRadius{inner_radius}};
            }
        };

        template <typename Ring_t>
        struct FlowRing : public Ring_t
        {
            using Ring_t::Ring_t;
            FlowRing(const Ring_t &r) : Ring_t{r} {}
        };
        template <typename Ring_t>
        struct TubeRing : public Ring_t
        {
            using Ring_t::Ring_t;
            TubeRing(const Ring_t &r) : Ring_t{r} {}
        };
        template <typename Ring_t>
        struct AnnulusRing : public Ring_t
        {
            using Ring_t::Ring_t;
            AnnulusRing(const Ring_t &r) : Ring_t{r} {}
        };
        template <typename Ring_t>
        struct ColumnRing : public Ring_t
        {
            using Ring_t::Ring_t;
            ColumnRing(const Ring_t &r) : Ring_t{r} {}
        };
        template <typename Ring_t>
        struct CementRing : public Ring_t
        {
            using Ring_t::Ring_t;
            CementRing(const Ring_t &r) : Ring_t{r} {}
        };

        bool is_tight_casing(
            const RealType lhs,
            const RealType rhs)
        {
            return std::abs(lhs - rhs) < 1e-12;
        }
        bool is_tight_casing(
            const Eigen::ArrayX<RealType> &lhs,
            const Eigen::ArrayX<RealType> &rhs)
        {
            assert(lhs.size() == rhs.size());

            const Eigen::ArrayX<RealType> temp{(lhs - rhs).abs()};
            return std::all_of(temp.cbegin(), temp.cend(), [](RealType v)
                               { return v < 1e-12; });
        }

        bool is_positive(const RealType v)
        {
            return v > 0.0;
        }
        bool is_positive(const Eigen::ArrayX<RealType> &v)
        {
            return std::all_of(v.cbegin(), v.cend(), [](RealType v)
                               { return v > 0.0; });
        }

        bool is_nan(const RealType v)
        {
            return std::isnan(v);
        }
        bool is_nan(const Eigen::ArrayX<RealType> &v)
        {
            return std::any_of(v.cbegin(), v.cend(), [](const RealType v)
                               { return std::isnan(v); });
        }

        template <typename Ring_t>
        struct Casing
        {
            using value_type = Ring_t::value_type;
            Casing(const std::vector<Ring_t> &completion)
                : sandwich{completion},
                  sandface_radius{completion.back().outer_radius},
                  flow_radius{completion.front().outer_radius},
                  column_outer_radius{completion[MaterialType::Column].outer_radius}
            {
                for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Cement; ++i)
                    assert(
                        is_tight_casing(
                            completion[i].inner_radius,
                            completion[i - 1ll].outer_radius));
            }

            const auto &back() const { return sandwich.back(); }
            const auto &front() const { return sandwich.front(); }

            const auto &flow() const { return front(); }

            const auto &operator[](auto i) const
            {
                return sandwich[i];
            }

            const value_type casing_volumetric_heat_capacity() const
            {
                // exclude "flow" at "i = 0" from summation!
                const auto start{MaterialType::Tube};
                value_type C{sandwich[start].linear_heat_capacity};
                for (ptrdiff_t i{start + 1}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    C += m.linear_heat_capacity;
                }
                C /= casing_area();
                return C;
            }

            const value_type cement_volumetric_heat_capacity() const
            {
                return sandwich[MaterialType::Cement].linear_heat_capacity / cement_area();
            }

            const value_type integral_casing_radial_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                const auto start{MaterialType::Tube};
                value_type L{sandwich[start].radial_heat_conductivity};
                for (ptrdiff_t i{start + 1}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    L += 1.0 / m.radial_heat_conductivity;
                }
                assert(is_positive(L));
                assert(!is_nan(L));
                return 1.0 / L;
            }

            const value_type integral_vertical_casing_heat_conductivity() const
            {
                // exclude "flow" at "i = 0"
                // as well as "cement2" at "i = end-1"
                // from summation!
                const auto start{MaterialType::Tube};
                value_type L{sandwich[start].integral_vertical_heat_conductivity};
                for (ptrdiff_t i{start + 1}; i <= MaterialType::Column; ++i)
                {
                    const auto &m = sandwich[i];
                    L += m.integral_vertical_heat_conductivity;
                }
                assert(is_positive(L));
                assert(!is_nan(L));
                return L / casing_area();
            }

            const value_type integral_vertical_cement_heat_conductivity() const
            {
                return sandwich[MaterialType::Cement].integral_vertical_heat_conductivity / cement_area();
            }

            const std::vector<Ring_t> sandwich;
            const value_type sandface_radius, column_outer_radius, flow_radius //, thickness
                ;

            const value_type area() const
            {
                const value_type out{cement_area() + casing_area()};
                return out;
            }

            const value_type cement_area() const
            {
                const value_type out{std::numbers::pi *
                                     (sandwich[MaterialType::Cement].thickness) *
                                     (sandwich[MaterialType::Cement].outer_radius + sandwich[MaterialType::Cement].inner_radius)};
                assert(is_positive(out));
                return out;
            }

            const value_type casing_area() const
            {
                const value_type out{std::numbers::pi *
                                     (sandwich[MaterialType::Column].outer_radius - sandwich[MaterialType::Tube].inner_radius) *
                                     (sandwich[MaterialType::Column].outer_radius + sandwich[MaterialType::Tube].inner_radius)};
                assert(is_positive(out));
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
            template <typename Container_t>
            ExtrudedRing(
                const StationaryPhaseProperties &props,
                const Container_t &inner_radius,
                const Container_t &thickness)
                : props{props},
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

            const StationaryPhaseProperties props;

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

        struct VarExtrudedRing
        {
            using value_type = Eigen::ArrayX<RealType>;

            VarExtrudedRing(
                const value_type &density,
                const value_type &specific_heat_capacity,
                const value_type &heat_conductivity,
                const value_type &inner_radius,
                const value_type &thickness)
                : thickness{thickness},
                  inner_radius{inner_radius},
                  outer_radius{inner_radius + thickness},
                  linear_heat_capacity{
                      linear_heat_capacity_calc(
                          inner_radius, thickness, density * specific_heat_capacity)},
                  radial_heat_conductivity{
                      radial_heat_conductivity_calc(
                          inner_radius, thickness, heat_conductivity)},
                  integral_vertical_heat_conductivity{
                      heat_conductivity * std::numbers::pi *
                      thickness * (thickness + 2.0 * inner_radius)}
            {
                for (auto id{0ll}; id < thickness.size(); ++id)
                    assert(std::abs(outer_radius(id) - inner_radius(id) - thickness(id)) < 1e-12);

                //        assert(thickness.size() == depth.size());
                assert(thickness.size() == inner_radius.size());
                assert(thickness.size() == outer_radius.size());
                assert(thickness.size() == density.size());
                assert(thickness.size() == specific_heat_capacity.size());
                assert(thickness.size() == heat_conductivity.size());
                assert(thickness.size() == linear_heat_capacity.size());
                assert(thickness.size() == radial_heat_conductivity.size());
                assert(thickness.size() == integral_vertical_heat_conductivity.size());
            }

            const Eigen::ArrayX<RealType>
                thickness, // depth,
                inner_radius, outer_radius,
                viscosity, specific_heat_capacity, heat_conductivity;

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
            static auto
            linear_heat_capacity_calc(
                const auto &inner_radius,
                const auto &thickness,
                const auto &volumetric_heat_capacity)
            {
                return std::numbers::pi *
                       thickness * (2 * inner_radius + thickness) *
                       volumetric_heat_capacity;
            }

            static auto
            radial_heat_conductivity_calc(
                const auto &inner_radius,
                const auto &thickness,
                const auto &heat_conductivity)
            {
                const auto temp{(1.0 + thickness / inner_radius).log()};
                Eigen::ArrayX<RealType> out{heat_conductivity / temp};
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