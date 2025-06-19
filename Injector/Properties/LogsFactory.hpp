#pragma once

#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/Factory.hpp>

namespace GPN
{
    namespace Logs
    {
        struct IsPermeableFactory
        {
            template <typename Grid_t>
            IsPermeableFactory(const Grid_t &grid)
                : IsPermeableFactory{
                      Logs::RawDataFactory::generate_is_permeable(
                          grid.dual_stencils),
                      grid}
            {
            }

            template <typename Container_t, typename Grid_t>
            IsPermeableFactory(const Container_t &is_permeable, const Grid_t &grid)
                : is_permeable{
                      IsPermeable{
                          StepPropertyGrid{
                              StepProperty{
                                  is_permeable},
                              grid}}}
            {
            }

            template <typename Container_t, typename Grid_t>
            static IsPermeable create(const Container_t &is_permeable, const Grid_t &grid)
            {
                return IsPermeable{
                    StepPropertyGrid{
                        StepProperty{
                            is_permeable},
                        grid}};
            }

            const auto &is_permeable_stencils() const
            {
                return is_permeable.log_vals;
            }

            IsPermeable is_permeable;
        };

        struct IsGhostLayerFactory
        {
            static IsGhostLayer create(
                const Logs::IsPermeable &is_permeable,
                const Logs::IsPerforated &is_perforated)
            {
                assert(is_perforated.size() == is_permeable.size());
                for (auto i{0ll}; i < (ptrdiff_t)is_perforated.size(); ++i)
                {
                    assert(
                        (is_perforated(i) == 0.0) ||
                        (is_perforated(i) == 1.0));
                    assert(
                        (is_perforated(i) == 0.0) ||
                        ((is_perforated(i) == 1.0) && (is_permeable(i) == 1.0)));
                }

                size_t predicate = 0ull;
                for (auto i{0ll}; i < (ptrdiff_t)is_perforated.size(); ++i)
                {
                    if ((is_perforated(i) == 0.0) &&
                        (is_permeable(i) == 1.0))
                    {
                        ++predicate;
                    }
                }
                // assume one or none ghost layers
                assert(predicate <= 1ull);

                return IsGhostLayer{
                    StepPropertyGrid{
                        StepPropertyContainer{
                            (is_permeable.log_vals - is_perforated.log_vals).eval()},
                        is_permeable.grid}};
            }

        };

        struct IsPerforatedFactory
        {
            template <typename Grid_t>
            static IsPerforated create(
                const auto &is_perforated,
                const auto &is_permeable,
                const Grid_t &grid)
            {
                assert(is_perforated.size() == is_permeable.size());
                for (auto i{0ll}; i < (ptrdiff_t)is_perforated.size(); ++i)
                {
                    assert(
                        (is_perforated[i] == 0.0) ||
                        (is_perforated[i] == 1.0));
                    assert(
                        (is_perforated[i] == 0.0) ||
                        ((is_perforated[i] == 1.0) && (is_permeable[i] == 1.0)));
                }

                size_t predicate = 0ull;
                for (auto i{0ll}; i < (ptrdiff_t)is_perforated.size(); ++i)
                {
                    if ((is_perforated[i] == 0.0) &&
                        (is_permeable[i] == 1.0))
                    {
                        ++predicate;
                    }
                }
                assert(predicate <= 1ull);

                // at least one perforated layer must exist
                assert(std::any_of(is_perforated.cbegin(), is_perforated.cend(), [](const RealType v){return v == 1.0;}));

                return IsPerforated{
                    StepPropertyGrid{
                        StepProperty{
                            is_perforated},
                        grid}};
            }

            const auto &is_permeable_stencils() const
            {
                return is_permeable.log_vals;
            }

            IsPerforated is_permeable;
        };

        struct PermeabilityFactory
        {
            template <typename Grid_t>
            static Permeability create(
                const auto &permeability,
                const auto &is_permeable,
                const Grid_t &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            permeability},
                        grid},
                    IsPermeableFactory::create(is_permeable, grid)};
            }
        };

        struct PorosityFactory
        {
            static Porosity create(
                const auto &porosity,
                const auto &is_permeable,
                const auto &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            porosity},
                        grid},
                    IsPermeableFactory::create(is_permeable, grid)};
            }
        };

        struct GeothermaFactory
        {
            /// @brief Create geotherma based on values table interpolation
            /// @param nodes Reference z-nodes for geotherma table
            /// @param vals Reference t-values for geotherma
            /// @param z_top Coordinate of the top
            /// @param q_grid Mesh nodes for temperature calculation
            /// @return
            static Geotherma create(
                const auto &nodes,
                const auto &vals,
                const RealType z_top,
                const auto &q_grid)
            {
                return {
                    StepPropertyGrid{
                        StepPropertyContainer{
                            interpolate(
                                nodes, vals, z_top, q_grid.mesh_nodes)},
                        q_grid}};
            }

            /// @brief Create const-value geotherms
            /// @param val Const temperature value
            /// @param q_grid Mesh nodes for temperature calculation
            /// @return
            static Geotherma create(
                const RealType val,
                const auto &q_grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            std::vector<RealType>(q_grid.mesh_nodes.size(), val)},
                        q_grid}};
            }

        private:
            static RealType interpolate_node(
                const auto &nodes, const auto &vals, RealType q_node, ptrdiff_t &left)
            {
                assert(left >= 0ll);
                // iterator to the point which is above the "node" value
                const auto it = std::upper_bound(nodes.cbegin() + left, nodes.cend(), q_node);
                const ptrdiff_t dist = std::distance(nodes.cbegin(), it);
                left = dist - 1ll;
                const ptrdiff_t right{dist};

                return (vals[right] - vals[left]) / (nodes[right] - nodes[left]) * (q_node - nodes[left]) + vals[left];
            }
            static auto interpolate(
                const auto &nodes, const auto &vals, const RealType z_top, const auto &q_nodes)
            {
                StepPropertyContainer out(q_nodes.size());

                ptrdiff_t left{0ll};
                for (auto i{0ll}; i < q_nodes.size(); ++i)
                {
                    out(i) = interpolate_node(
                        nodes, vals, z_top + q_nodes[i], left);
                }
                return out;
            }
        };

        struct SkinFactory
        {
            static SkinFactor create(
                const auto &skin,
                const auto &is_permeable,
                const auto &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            skin},
                        grid},
                    IsPermeableFactory::create(is_permeable, grid)};
            }
        };

        struct ExtPressureFactory
        {
            static ExternalPressure create(
                const auto &pressure,
                const auto &is_permeable,
                const auto &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            pressure},
                        grid},
                    IsPermeableFactory::create(is_permeable, grid)};
            }
        };

        struct HeatConductivityFactory
        {
            static HeatConductivity create(
                const auto &heat_conductivity,
                const auto &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            heat_conductivity},
                        grid}};
            }
        };

        struct SolidDensityFactory
        {
            static SolidDensity create(
                const auto &solid_density,
                const auto &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            solid_density},
                        grid}};
            }
        };

        struct ThermalDiffusivityFactory
        {
            static ThermalDiffusivity create(
                const auto &capacity,
                const auto &conductivity,
                const auto &grid)
            {
                return {StepPropertyGrid{
                            StepProperty{conductivity / capacity}},
                        grid};
            }
        };

        struct SolidSpecificHeatCapacityFactory
        {
            static SolidSpecificHeatCapacity create(
                const auto &solid_specific_heatcapacity,
                const auto &grid)
            {
                assert(solid_specific_heatcapacity.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                return {
                    StepPropertyGrid{
                        StepProperty{
                            solid_specific_heatcapacity},
                        grid}};
            }
        };

        struct SolidVolumetricHeatCapacityFactory
        {
            static SolidVolumetricHeatCapacity create(
                const auto &solid_volumetric_heatcapacity,
                const auto &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            solid_volumetric_heatcapacity},
                        grid}};
            }

            static SolidVolumetricHeatCapacity create(
                const auto &density,
                const auto &heat_capacity,
                const auto &grid)
            {
                return {
                    SolidDensityFactory::create(density, grid),
                    SolidSpecificHeatCapacityFactory::create(heat_capacity, grid)};
            }
        };

        struct MediumHeatVolumetricCapacityFactory
        {
            static MediumHeatVolumetricCapacity create(
                const auto &porosity,
                const auto &solid_vol_heatcapacity,
                const auto &fluid,
                const auto &grid)
            {
                return {StepPropertyGrid{StepProperty{
                                             porosity * fluid.volumetric_heat_capacity +
                                             (1.0 - porosity) * solid_vol_heatcapacity},
                                         grid}};
            }

            static MediumHeatVolumetricCapacity create(
                const auto &porosity,
                const auto &solid_density,
                const auto &solid_heat_capacity,
                const auto &fluid,
                const auto &grid)
            {
                return create(porosity, solid_density * solid_heat_capacity, fluid, grid);
            }
        };

        struct RFPFactory
        {
            template <typename Container_t, typename IsPermeable_t>
            static auto create_from_container(
                const Container_t &rfp,
                const IsPermeable_t &is_permeable)
            {
                return RFP{
                    StepPropertyGrid{rfp, is_permeable.grid},
                    is_permeable};
            }

            template <typename Record_t, typename Well_t>
            static auto create_from_well(
                const Record_t& history_record,
                const Well_t &well)
            {
                return RFP{
                    StepPropertyGrid{well.get_RFP(history_record), well.is_permeable.grid},
                    well.is_permeable};
            }
        };

        struct WFPFactory
        {
            template <typename Container_t, typename IsPerforated_t>
            static auto create_from_container(
                const Container_t &wfp,
                const IsPerforated_t &is_perforated)
            {
                return WFP{
                    StepPropertyGrid{wfp, is_perforated.grid},
                    is_perforated};
            }

            template <typename Record_t, typename Well_t>
            static auto create_from_well(
                const Record_t history_record,
                const Well_t &well)
            {
                return WFP{
                    StepPropertyGrid{well.get_WFP(history_record), well.is_perforated.grid},
                    well.is_perforated};
            }
        };

        struct HydrodynamicLogsFactory : public IsPermeableFactory
        {
            template <typename Grid_t>
            HydrodynamicLogsFactory(const Grid_t &grid)
                : HydrodynamicLogsFactory{IsPermeableFactory{grid}, grid}
            {
            }

            template <typename Grid_t>
            HydrodynamicLogsFactory(
                const IsPermeableFactory &is_permeable_factory,
                const Grid_t &grid)
                : HydrodynamicLogsFactory{
                      is_permeable_factory,
                      RawDataFactory::generate_porosity(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      RawDataFactory::generate_permeability(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      RawDataFactory::generate_ext_pressure(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      RawDataFactory::generate_skin(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      grid}
            {
            }

            HydrodynamicLogsFactory(
                const IsPermeableFactory &is_permeable_factory,
                const auto &porosity,
                const auto &permeability,
                const auto &ext_pressure,
                const auto &skin,
                const auto &grid)
                : IsPermeableFactory{is_permeable_factory},
                  permeability{
                      StepPropertyGrid{
                          StepProperty{permeability},
                          grid},
                      is_permeable_factory.is_permeable},
                  porosity{
                      StepPropertyGrid{
                          StepProperty{porosity},
                          grid},
                      is_permeable_factory.is_permeable},
                  ext_pressure{
                      StepPropertyGrid{
                          StepProperty{ext_pressure},
                          grid},
                      is_permeable_factory.is_permeable},
                  skin{
                      StepPropertyGrid{
                          StepProperty{skin},
                          grid},
                      is_permeable_factory.is_permeable}
            {
            }

            template <typename Container_t, typename Grid_t>
            HydrodynamicLogsFactory(
                const Container_t &is_permeable,
                const Container_t &porosity,
                const Container_t &permeability,
                const Container_t &ext_pressure,
                const Container_t &skin,
                const Grid_t &grid)
                : HydrodynamicLogsFactory{
                      IsPermeableFactory{is_permeable, grid},
                      porosity, permeability, ext_pressure, skin,
                      grid}
            {
            }

            const Permeability permeability;
            const Porosity porosity;
            const ExternalPressure ext_pressure;
            const SkinFactor skin;
        };

        struct HeatLogsFactory
        {
            template <typename Container_t, typename Grid_t>
            HeatLogsFactory(
                const Container_t &conductivity,
                const Grid_t &grid)
                : conductivity{
                      Logs::StepPropertyGrid{
                          Logs::StepProperty{
                              conductivity},
                          grid}}
            {
            }

            const HeatConductivity conductivity;
        };
    } // Logs
} // GPN