#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>

namespace GPN
{
    namespace DataContainers
    {
        struct SpatialGridContainer
        {
            using VR = std::vector<RealType>;
            SpatialGridContainer(
                const VR &x1_stencils,
                const VR &x2_stencils)
            {
            }

        protected:
            const VR x1_stencils;
            const VR x2_stencils;
        };

        struct CylinderGridContainer : public SpatialGridContainer
        {
            using SpatialGridContainer::SpatialGridContainer;

            const auto& z_stencils() const{
                return x1_stencils;
            }
            const auto& r_stencils() const{
                return x2_stencils;
            }
        }

    } // DataContainers

} // GPN