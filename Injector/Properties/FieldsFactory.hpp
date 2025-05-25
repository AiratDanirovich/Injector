#pragma once

// #include <Injector/Grids/Grids2D.hpp>
#include <Injector/Properties/PhysicalField.hpp>
#include <Injector/Properties/LogsFactory.hpp>

namespace GPN
{
    namespace Properties
    {
        struct Factory
        {
            // static auto generate_heatconductivity_Property(
            //     const auto &dual_stencils,
            //     const auto grid)
            // {
            //     return Properties::HeatConductivity{
            //         Logs::HeatConductivity{
            //             Logs::StepPropertyGrid{
            //                 Logs::StepProperty{dual_stencils},
            //                 grid->first_coord}},
            //         grid};
            // }

            static auto generate_solid_volumetric_heatcapacity_Property(
                const auto &soliddensity_stencils,
                const auto &solid_specific_heatcapacity_stencils,
                const auto grid)
            {
                const auto solid_density{
                    Logs::SolidDensityFactory::create(
                        soliddensity_stencils, grid->first_coord)};

                const auto solid_specific_heatcapacity{
                    Logs::SolidSpecificHeatCapacityFactory::create(
                        solid_specific_heatcapacity_stencils, grid->first_coord)};

                return Logs::SolidVolumetricHeatCapacity{
                    solid_density, solid_specific_heatcapacity};
            }

            // static auto generate_volumetric_heatcapacity_Property(
            //     const auto &is_permeable_stencils,
            //     const auto &porosity_stencils,
            //     const auto &soliddensity_stencils,
            //     const auto &solid_specific_heatcapacity_stencils,
            //     const auto &fluid,
            //     const auto grid)
            // {
            //     const auto &z_grid{grid->first_coord};
            //     const auto is_permeable{
            //         Logs::IsPermeable{
            //             Logs::StepPropertyGrid{
            //                 Logs::StepProperty{
            //                     is_permeable_stencils},
            //                 z_grid}}};

            //     const auto porosity{
            //         Logs::Porosity{
            //             Logs::StepPropertyGrid{
            //                 Logs::StepProperty{
            //                     porosity_stencils},
            //                 z_grid} * // guarantee that porosity is zero in rocks
            //                 is_permeable,
            //             is_permeable}};

            //     const auto solid_volumetric_heatcapacity{
            //         Properties::Factory::generate_solid_volumetric_heatcapacity_Property(
            //             soliddensity_stencils,
            //             solid_specific_heatcapacity_stencils,
            //             grid)};

            //     const auto capacity{
            //         Logs::MediumHeatVolumetricCapacityFactory::create(
            //             porosity, solid_volumetric_heatcapacity, fluid, grid->first_coord)};

            //     return Properties::MediumHeatVolumetricCapacity{
            //         capacity, grid};
            // }
        };

    } // Properties
} // GPN