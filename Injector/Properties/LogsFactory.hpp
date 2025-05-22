#pragma once

#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    namespace Logs
    {
        struct IsPermeableFactory
        {
            template <typename Grid_t>
            IsPermeableFactory(const Grid_t &grid)
                : IsPermeableFactory{
                      Logs::Factory::generate_is_permeable_StepProperty(
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

            const auto &is_permeable_stencils() const
            {
                return is_permeable.log_vals;
            }

            IsPermeable is_permeable;
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
                const IsPermeableFactory& is_permeable_factory, 
                const Grid_t &grid)
                : HydrodynamicLogsFactory{is_permeable_factory,
                      Factory::generate_porosity_StepProperty(
                          grid.dual_stencils,
                          is_permeable_factory.is_permeable_stencils()),
                      Factory::generate_permeability_StepProperty(
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
                      is_permeable_factory.is_permeable}
            {
            }

            template <typename Container_t, typename Grid_t>
            HydrodynamicLogsFactory(
                const Container_t &is_permeable,
                const Container_t &porosity,
                const Container_t &permeability,
                const Grid_t &grid)
                : HydrodynamicLogsFactory{
                    IsPermeableFactory{is_permeable, grid}, 
                    porosity, permeability, 
                    grid}
            {
            }

            const Permeability permeability;
            const Porosity porosity;
        };

    } // Logs

} // GPN