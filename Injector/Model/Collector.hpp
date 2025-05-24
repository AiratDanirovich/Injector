#pragma once

#include <Injector/Properties/LogsFactory.hpp>

namespace GPN
{
    namespace Logs
    {
        namespace Rocks
        {
            struct CoreSampleData
            {
                template <typename Container_t, typename Grid_t>
                CoreSampleData(
                    const Container_t &is_permeable_stencils,
                    const Container_t &porosity_stencils,
                    const Container_t &permeability_stencils,
                    const Grid_t &grid)
                    : is_permeable{IsPermeableFactory::create(is_permeable_stencils, grid)},
                      permeability{PermeabilityFactory::create(permeability_stencils, is_permeable_stencils, grid)},
                      porosity{PorosityFactory::create(permeability_stencils, is_permeable_stencils, grid)}
                {
                }

                IsPermeable is_permeable;
                Permeability permeability;
                Porosity porosity;
            };
        } // Rocks

        namespace Hyrodynamics
        {
            struct Hyrodynamics
            {
            };

        } // Hydrodynamics

    } // Properties
} // GPN