#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids.h>

namespace GPN
{
    namespace Grids
    {
        struct Factory
        {
            static auto generate_stencils_uniform(RealType a, RealType b, ptrdiff_t n)
            {
                RealType step{(b - a) / (n - 1)};

                std::vector<RealType> result(n);

                result[0] = a;
                for (ptrdiff_t idx{1}; idx < n; ++idx)
                    result[idx] = result[idx - 1] + step;

                return result;
            }

            static auto create_grid_2D(ptrdiff_t n)
            {
                auto stencils{Factory::generate_stencils_uniform(0, 1, n)};

                auto nodes{GridNodes{stencils}};

                auto z_grid{
                    Grid1D<CoordinateTypes::Z>{nodes}};

                auto r_grid{
                    Grid1D<CoordinateTypes::RadialCylinderCoordinate>{nodes}};

                return 
                    StructuredCylinderGrid2D{
                        Grid1D<CoordinateTypes::Z>{nodes},
                        Grid1D<CoordinateTypes::RadialCylinderCoordinate>{nodes}};
            }
        };
    }
}