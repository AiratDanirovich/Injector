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
                : is_permeable{
                      IsPermeable{
                          StepPropertyGrid{
                              StepProperty{
                                  Logs::Factory::generate_is_permeable_StepProperty(grid.dual_stencils)},
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
                : IsPermeableFactory{grid},
                  permeability{
                      Permeability{
                          StepPropertyGrid{
                              StepProperty{
                                  Factory::generate_permeability_StepProperty(
                                      grid,
                                      is_permeable_stencils())},
                              grid},
                          is_permeable}},
                  porosity{
                      Porosity{
                          StepPropertyGrid{
                              StepProperty{
                                  Factory::generate_porosity_StepProperty(
                                      grid,
                                      is_permeable_stencils())},
                              grid},
                          is_permeable}}
            {
            }
            const Permeability permeability;
            const Porosity porosity;
        };

    } // Logs

} // GPN