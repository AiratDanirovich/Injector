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

        struct PermeabilityFactory
        {
            template <typename Container_t, typename Grid_t>
            static Permeability create(
                const Container_t &permeability,
                const Container_t &is_permeable,
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
                return {StepPropertyGrid{
                    porosity * fluid.volumetric_heat_capacity +
                        (1.0 - porosity) * solid_vol_heatcapacity,
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
            static auto create(
                const Container_t &rfp,
                const IsPermeable_t &is_permeable)
            {
                return RFP{
                    StepPropertyGrid{rfp, is_permeable.grid},
                    is_permeable};
            }

            template <typename Well_t>
            static auto create(
                RealType well_rate,
                const Well_t &well)
            {
                return RFP{
                    StepPropertyGrid{well.get_RFP(well_rate), well.is_permeable.grid},
                    well.is_permeable};
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
                      RawDataFactory::generate_skin_StepProperty(
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