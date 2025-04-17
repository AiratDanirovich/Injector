#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/SpacialGrids.hpp>
#include <Injector/Grids/Grids2D.hpp>

namespace GPN
{
    namespace Grids
    {
        struct Factory
        {
            static auto generate_dual_grid_stencils_uniform(RealType a, RealType b, ptrdiff_t n)
            {
                assert(n > 1);
                RealType step{(b - a) / (n - 1)};

                std::vector<RealType> result(n);

                result.front() = a;
                for (ptrdiff_t idx{1}; idx < n-1; ++idx)
                    result[idx] = result[idx - 1] + step;
                result.back() = b;

                return result;
            }

            static auto create_cartesian_grid_2D(ptrdiff_t n)
            {
                auto stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, n)};

                auto nodes{SpacialGrid{stencils}};

                auto x_grid{
                    AxesGrid<CoordinateTypes::X>{nodes}};

                auto y_grid{
                    AxesGrid<CoordinateTypes::Y>{nodes}};

                return 
                    StructuredXYGrid2D{x_grid, y_grid};
            }
            
            static auto create_cylinder_grid_2D(ptrdiff_t n)
            {
                auto stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, n)};

                auto nodes{GridDual{stencils}};

                auto z_grid{
                    AxesGrid<CoordinateTypes::Z>{nodes}};

                auto r_grid{
                    AxesGrid<CoordinateTypes::R_CylCoord>{nodes}};

                return 
                    StructuredCylinderGrid2D{z_grid,r_grid};
            }
        };
    }
}