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
                      Logs::RawdataFactory::generate_is_permeable_StepProperty(
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
            template <typename Container_t, typename Grid_t>
            static Porosity create(
                const Container_t &porosity,
                const Container_t &is_permeable,
                const Grid_t &grid)
            {
                return {
                    StepPropertyGrid{
                        StepProperty{
                            porosity},
                        grid},
                    IsPermeableFactory::create(is_permeable, grid)};
            }
        };

        struct RFPFactory
        {
            template <typename Container_t, typename IsPermeable_t>
            static RFP create(
                const Container_t &rfp,
                const IsPermeable_t &is_permeable)
            {
                return RFP{
                    StepPropertyGrid{rfp, is_permeable.grid},
                    is_permeable};
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
                      RawdataFactory::generate_porosity_StepProperty(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      RawdataFactory::generate_permeability_StepProperty(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      RawdataFactory::generate_ext_pressure_StepProperty(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      RawdataFactory::generate_skin_StepProperty(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      grid}
            {
            }

            template <typename Container_t, typename Grid_t>
            HydrodynamicLogsFactory(
                const IsPermeableFactory &is_permeable_factory,
                const Container_t &porosity,
                const Container_t &permeability,
                const Container_t &ext_pressure,
                const Container_t &skin,
                const Grid_t &grid)
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