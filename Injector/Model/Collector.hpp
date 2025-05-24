#pragma once

#include <Injector/Properties/LogsFactory.hpp>

namespace GPN
{
    namespace Logs
    {
        namespace Rocks
        {
            struct CoreSampleLogs
            {
                CoreSampleLogs(
                    const auto &is_permeable_stencils,
                    const auto &porosity_stencils,
                    const auto &permeability_stencils,
                    const auto &grid)
                    : is_permeable{IsPermeableFactory::create(is_permeable_stencils, grid)},
                      permeability{PermeabilityFactory::create(permeability_stencils, is_permeable_stencils, grid)},
                      porosity{PorosityFactory::create(permeability_stencils, is_permeable_stencils, grid)}
                {
                }

                IsPermeable is_permeable;
                Permeability permeability;
                Porosity porosity;
            };

            struct HeatLogs
            {
                HeatLogs(
                    const auto &solid_density,
                    const auto &solid_specific_heatcapacity,
                    const auto &heat_conductivity,
                    const auto &grid)
                    : solid_density{
                        SolidDensityFactory::create(solid_density, grid)},
                      solid_specific_heatcapacity{
                        SolidSpecificHeatCapacityFactory::create(solid_specific_heatcapacity, grid)},
                      heat_conductivity{
                        HeatConductivityFactory::create(heat_conductivity, grid)},
                      solid_vol_heatcapacity{
                        SolidVolumetricHeatCapacityFactory::create(solid_density, solid_specific_heatcapacity, grid)}
                {
                }

                const SolidDensity solid_density;
                const SolidSpecificHeatCapacity solid_specific_heatcapacity;
                const HeatConductivity heat_conductivity;
                const SolidVolumetricHeatCapacity solid_vol_heatcapacity;

            private:
                template <typename T>
                static T multiply(
                    const T &lhs,
                    const T &rhs)
                {
                    T out(lhs.size(), 0.0);
                    for (auto i{0ll}; i < lhs.size(); ++i)
                        out[i] = lhs[i] * rhs[i];

                    return out;
                }
            };
        } // Rocks

        namespace Hydrodynamics
        {
            struct Hydrodynamics
            {
                Hydrodynamics(
                    const auto &is_permeable_stencils,
                    const auto &ext_pressure,
                    const auto &skin,
                    const auto &grid)
                    : skin{SkinFactory::create(skin, is_permeable_stencils, grid)},
                      ext_pressure{ExtPressureFactory::create(ext_pressure, is_permeable_stencils, grid)}
                {
                }

                const ExternalPressure ext_pressure;
                const SkinFactor skin;
            };

        } // Hydrodynamics
    } // Properties
} // GPN