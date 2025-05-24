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
                    const auto &porosity_stencils,
                    const auto &permeability_stencils,
                    const auto &grid)
                    : solid_density{},
                      solid_specific_heatcapacity{},
                      solid_vol_heatcapacity{}
                {
                }

                const SolidDensity solid_density;
                const SolidSpecificHeatCapacity solid_specific_heatcapacity;
                const SolidVolumetricHeatCapacity solid_vol_heatcapacity;
                const HeatConductivity heat_conductivity;
            };
        } // Rocks

        namespace Hyrodynamics
        {
            struct Hyrodynamics
            {
                Hyrodynamics(
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