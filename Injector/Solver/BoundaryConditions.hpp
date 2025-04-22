#pragma once

#include <memory>
#include <Eigen/Core>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids2D.hpp>

namespace GPN
{
    namespace BoundaryConditions
    {
        struct BCSouthNorth
        {
            BCSouthNorth(
                const BCSouth &south,
                const BCNorth &north,
                const Grids::GridDual &grid,
                std::shared_ptr<const BCFunctorBase> functor,
                RealType t = 0.0)
                : south{south}, north{north},
                  grid{grid},
                  south_vals(grid.mesh_size()),
                  north_vals(grid.mesh_size()),
                  functor{functor}
            {
                assert(south_vals.size() == north_vals.size());
                set_vals(t);
            }

            BCSouth south;
            BCNorth north;
            const Grids::GridDual &grid;

            Eigen::ArrayX<RealType> south_vals, north_vals;
            std::shared_ptr<const BCFunctorBase> functor;

            void set_vals(RealType t)
            {
                for (std::ptrdiff_t i{0ull}; i < south_vals.size(); ++i)
                {
                    south_vals[i] =
                        (*functor)(south.fixed_x, grid.coordinate(i), t);
                    north_vals[i] =
                        (*functor)(north.fixed_x, grid.coordinate(i), t);
                }
            }
        };

        struct BCEastWest
        {
            BCEastWest(const BCEast &east,
                       const BCWest &west,
                       const Grids::GridDual &grid,
                       std::shared_ptr<const BCFunctorBase> functor,
                       RealType t = 0.0)
                : east{east},
                  west{west},
                  grid{grid},
                  east_vals(grid.mesh_size()),
                  west_vals(grid.mesh_size()),
                  functor{functor}
            {
                assert(east_vals.size() == west_vals.size());
                set_vals(t);
            }

            BCEast east;
            BCWest west;
            const Grids::GridDual &grid;

            // f(x,y) for the first-type boundary condition
            std::shared_ptr<const BCFunctorBase> functor;

            Eigen::ArrayX<RealType> east_vals, west_vals;

            // set values u(x,y) at fixed y = y_east and y = y_west
            void set_vals(RealType t)
            {
                for (std::ptrdiff_t i{0ull}; i < east_vals.size(); ++i)
                {
                    east_vals[i] = (*functor)(grid.coordinate(i), east.fixed_y, t);
                    west_vals[i] = (*functor)(grid.coordinate(i), west.fixed_y, t);
                }
            }
        };

        struct BoundaryConditions
        {
            template <typename Grid_t>
            BoundaryConditions(
                const Grid_t &grid,
                std::shared_ptr<const BCFunctorBase> functor)
                : south_north{
                      BCSouth{grid.first_coord.front()},
                      BCNorth{grid.first_coord.back()},
                      grid.second_coord,
                      functor},
                  east_west{
                        BCEast{grid.second_coord.back()}, 
                        BCWest{grid.second_coord.front()}, 
                        grid.first_coord, 
                        functor
                    }
            {
            }

            BoundaryConditions(const BoundaryConditions &) = default;

            RealType east_vals(auto i) const
            {
                return east_west.east_vals[i];
            }
            RealType west_vals(auto i) const
            {
                return east_west.west_vals[i];
            }

            RealType south_vals(auto i) const
            {
                return south_north.south_vals[i];
            }
            RealType north_vals(auto i) const
            {
                return south_north.north_vals[i];
            }

            void set_vals(RealType t)
            {
                east_west.set_vals(t);
                south_north.set_vals(t);
            }

        protected:
            BCEastWest east_west;
            BCSouthNorth south_north;
        };
    } // BoundaryConditions
} // GPN
