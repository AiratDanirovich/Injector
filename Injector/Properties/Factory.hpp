#pragma once

#include <vector>

#include <Injector/Grids/Grids1D.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>

namespace GPN
{
    namespace Logs
    {
        struct Factory
        {
            using Grid_t = Grids::AxesGrid<CoordinateTypes::Z>;

            template <typename Container_t>
            static auto generate_permeability_StepProperty(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 500 : 300;
                return vals;
            }

            template <typename Container_t>
            static auto generate_porosity_StepProperty(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 0.2 : 0.5;
                return vals;
            }

            template <typename Container_t>
            static auto generate_solid_density_StepProperty(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 800.0 : 1000.0;
                return vals;
            }

            template <typename Container_t>
            static auto generate_solid_specific_heatcapacity_StepProperty(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 0.7 : 0.9;
                return vals;
            }

            template <typename Container_t>
            static auto generate_is_permeable_StepProperty(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 0.0 : 1.0;
                return vals;
            }

            template <typename Container_t>
            static auto generate_conductivity_StepProperty(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 200 : 1000;
                return vals;
            }

            static auto generate_rates_StepProperty(const Grids::GridDualStencils &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 50.0 : 0.0;
                return vals;
            }

            template <typename Container_t>
            static auto generate_soliddensity_Log(
                const Container_t &dual_stencils,
                const auto grid)
            {
                return Logs::SolidDensity{
                    Logs::StepPropertyGrid{
                        Logs::StepProperty{
                            dual_stencils},
                        grid->first_coord}};
            }

            template <typename Container_t>
            static auto generate_solid_specific_heatcapacity_Log(
                const Container_t &solid_specific_heatcapacity_stencils,
                const auto grid)
            {
                return Logs::SolidSpecificHeatCapacity{
                    StepPropertyGrid{
                        StepProperty{
                            solid_specific_heatcapacity_stencils},
                        grid->first_coord}};
            }
        };
    }

    namespace Properties
    {
        struct Factory
        {
            template <typename Container_t>
            static Properties::HeatConductivity generate_heatconductivity_Property(
                const Container_t &dual_stencils,
                const auto grid)
            {
                const auto t{Logs::StepProperty{dual_stencils}};
                const auto v{Logs::StepPropertyGrid{
                    t,
                    grid->first_coord}};
                const auto b{Logs::HeatConductivity{
                    v}};

                return Properties::HeatConductivity{b, grid};
            }

            template <typename Container_t>
            static auto generate_solid_volumetric_heatcapacity_Property(
                const Container_t &soliddensity_stencils,
                const Container_t &solid_specific_heatcapacity_stencils,
                const auto grid)
            {
                const auto solid_density{
                    Logs::Factory::generate_soliddensity_Log(
                        soliddensity_stencils, grid)};

                const auto solid_specific_heatcapacity{
                    Logs::Factory::generate_solid_specific_heatcapacity_Log(
                        solid_specific_heatcapacity_stencils, grid)};

                return Logs::SolidVolumetricHeatCapacity{
                    solid_density, solid_specific_heatcapacity};
            }

            template <typename Container_t>
            static auto generate_volumetric_heatcapacity_Property(
                const Container_t &is_permeable_stencils,
                const Container_t &porosity_stencils,
                const Container_t &soliddensity_stencils,
                const Container_t &solid_specific_heatcapacity_stencils,
                const GPN::PhaseProperties &fluid,
                const auto grid)
            {
                const auto &z_grid{grid->first_coord};
                const auto is_permeable{
                    Logs::IsPermeable{
                        Logs::StepPropertyGrid{
                            Logs::StepProperty{
                                is_permeable_stencils},
                            z_grid}}};

                const auto porosity{
                    Logs::Porosity{
                        Logs::StepPropertyGrid{
                            Logs::StepProperty{
                                porosity_stencils},
                            z_grid} * // guarantee that porosity is zero in rocks
                            is_permeable,
                        is_permeable}};

                const auto solid_volumetric_heatcapacity{
                    Properties::Factory::generate_solid_volumetric_heatcapacity_Property(
                        soliddensity_stencils,
                        solid_specific_heatcapacity_stencils,
                        grid)};

                const auto capacity{
                    Logs::HeatVolumetricCapacity{
                        porosity, solid_volumetric_heatcapacity, fluid}};

                return Properties::HeatVolumetricCapacity{
                    capacity, grid};
            }
        };

    } // Properties
}