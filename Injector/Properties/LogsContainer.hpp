#pragma once

#include <map>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    namespace Logs
    {
        enum struct Log
        {
            is_permeable,
            permeability,
            porosity,
            solid_density,
            solid_specific_heat_capacity,
            solid_volumetric_heat_capacity,
            heat_volumetric_capacity,
            ext_pressure
        };
        
        template <typename Grid_1D_t>
        struct GeologyLogsContainer
        {
            GeologyLogsContainer(
                const Grid_1D_t &grid,
                const PhaseProperties &fluid,
                const MeshNodesContainer &is_permeable,
                const MeshNodesContainer &permeability,
                const MeshNodesContainer &porosity,
                const MeshNodesContainer &solid_density,
                const MeshNodesContainer &solid_specific_heat_capacity,
                const MeshNodesContainer &ext_pressure)
                : grid{grid}
            {
                const SolidDensity solid_density_log{
                    StepPropertyGrid{
                        StepProperty{solid_density},
                        grid}};
                logs.insert(
                    {Log::solid_density,
                     solid_density_log.log_vals});

                const SolidSpecificHeatCapacity solid_specific_heat_capacity_log{
                    StepPropertyGrid{
                        StepProperty{solid_specific_heat_capacity},
                        grid}};
                logs.insert(
                    {Log::solid_specific_heat_capacity,
                     solid_specific_heat_capacity_log.log_vals});

                const SolidVolumetricHeatCapacity solid_volumetric_heat_capacity_log{
                    solid_density_log, solid_specific_heat_capacity_log};
                logs.insert(
                    {Log::solid_volumetric_heat_capacity,
                     solid_volumetric_heat_capacity_log.log_vals});

                const IsPermeable is_permeable_log{
                    IsPermeable{
                        StepPropertyGrid{
                            StepProperty{
                                is_permeable},
                            grid}}};

                logs.insert(
                    {Log::is_permeable, is_permeable_log.log_vals});

                logs.insert(
                    {Log::permeability,
                     Permeability{
                         StepPropertyGrid{
                             StepProperty{permeability},
                             grid},
                         is_permeable_log}
                         .log_vals});

                const Porosity porosity_log{
                    StepPropertyGrid{
                        StepProperty{porosity},
                        grid},
                    is_permeable_log};
                logs.insert({Log::porosity,
                             porosity_log.log_vals});

                HeatVolumetricCapacity heat_volumetric_capacity_log(
                    porosity_log,
                    solid_volumetric_heat_capacity_log,
                    fluid);

                logs.insert(
                    {Log::heat_volumetric_capacity,
                     ExternalPressure{
                         StepPropertyGrid{
                             StepProperty{ext_pressure},
                             grid},
                         is_permeable_log}
                         .log_vals});





                logs.insert(
                    {Log::ext_pressure,
                     ExternalPressure{
                         StepPropertyGrid{
                             StepProperty{ext_pressure},
                             grid},
                         is_permeable_log}
                         .log_vals});
            }

            std::map<Log, const MeshNodesContainer> logs;

        private:
            const Grid_1D_t &grid;
        };

    } // Logs

} // GPN