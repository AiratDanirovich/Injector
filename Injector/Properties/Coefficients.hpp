#pragma once

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
// #include <Injector/Grids/UniformGrid2D.h>

namespace GPN
{
    namespace EqSolver
    {
        namespace Properties
        {
            struct Fields
            {
                using Capacity_f = CellVolumeContainer2D;
                using Source_f = CellVolumeContainer2D;
                using Conductivity_f = CellVolumeContainer2D;

                Capacity_f capacity_vol_f;
                Source_f source_vol_f;
                Conductivity_f conductivity_f;

                // template <
                //     typename Grid_t,
                //     typename Factory_t
                // >
                // Fields(const Grid_t &grid,
                //        const Factory_t &factory)
                //     : Fields{grid,
                //              Factory_t::Capacity{},
                //              Factory_t::Conductivity{},
                //              Factory_t::Source{}
                //     }
                // {
                // }

                template <
                    typename Grid_t,
                    typename Capacity_t,
                    typename Conductivity_t,
                    typename Source_t
                >
                Fields(const Grid_t &grid,
                       const Capacity_t &capacity,
                       const Conductivity_t &conductivity,
                       const Source_t &source)
                    : capacity_vol_f{set_functor(grid, capacity) * grid.volumes()},
                      conductivity_f{set_functor(grid, conductivity)},
                      source_vol_f{set_functor(grid, source) * grid.volumes()}
                {
                }

            protected:
                template <
                    typename Grid_t, 
                    typename Functor
                >
                static auto set_functor(
                    const Grid_t &grid,
                    const Functor &functor)
                {
                    Eigen::ArrayXX<RealType> functor_nodes{
                        grid.first_coord.size(),
                        grid.second_coord.size()};

#pragma omp parallel for
                    for (std::ptrdiff_t j = 0; j < functor_nodes.cols(); ++j)
                        for (std::ptrdiff_t i = 0; i < functor_nodes.rows(); ++i)
                        {
                            functor_nodes(i, j) =
                                functor(
                                    grid.first_coord.coord(i),
                                    grid.second_coord.coord(j)
                                );
                        }

                    return functor_nodes;
                }
            };
        } // Properties
    } // EqSolver
} // GPN