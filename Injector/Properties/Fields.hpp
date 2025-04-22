#pragma once

#include <Injector/Properties/Logs.hpp>
#include <Injector/Grids/Grids2D.hpp>


namespace GPN
{
    namespace Field2D
    {



                    /**
                 * @brief Interpolate conductivity in y-direction, 
                 * to get values on y-boundaries of control volumes. 
                 * Each volume has two boundaries in a fixed direction, 
                 * but two boundaries at the domain boundary are excluded, 
                 * and the number of boundaries in container is
                 * by one less than the number of control volumes. 
                 * @tparam Grid_t
                 * @param conductivity_nodes
                 * @param grid
                 * @return auto
                 */
//                 template <typename Grid_t>
//                 auto set_conductivity_y_bounds(
//                     const Conductivity_f &conductivity_nodes,
//                     const Grid_t &grid) const
//                 {
//                     Conductivity_f conductivity_y_bounds{
//                         grid.first_coord.size(),
//                         grid.second_coord.size() - 1};

// #pragma omp parallel for
//                     for (ptrdiff_t i = 0; i < conductivity_y_bounds.rows(); ++i)
//                         for (ptrdiff_t j = 0; j < conductivity_y_bounds.cols(); ++j)
//                             conductivity_y_bounds(i, j) =
//                                 (conductivity_nodes(i, j + 1) + conductivity_nodes(i, j)) / 2.0;
//                     return conductivity_y_bounds;
//                 }
                
    } // Fields2D
} // GPN