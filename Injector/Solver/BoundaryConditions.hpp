#pragma once

#include <memory>
#include <Eigen/Core>
#include <Eigen/SparseCore>

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
                for (std::ptrdiff_t i{0ll}; i < south_vals.size(); ++i)
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

        struct MatrixView
        {
            using MatrixRow_t = Eigen::Block<Eigen::SparseMatrix<RealType>, 1, -1, false>;
            using RHS_t = Eigen::Block<Eigen::VectorXd, 1, 1, false>;
            MatrixView(MatrixRow_t A,
                       RHS_t rhs,
                       std::ptrdiff_t diag,
                       std::ptrdiff_t neib)
                : matrix{A}, rhs{rhs}, diag{diag}, neib{neib}
            {
            }

            MatrixRow_t matrix;
            RHS_t rhs;
            std::ptrdiff_t diag, neib;
        };

        struct BoundaryConditions
        {
            template <typename Grid_t>
            BoundaryConditions(
                const Grid_t &grid,
                std::shared_ptr<const BCFunctorBase> functor,
                BoundaryCondition::BCType bc_type = BoundaryCondition::first)
                : south_north{
                      BCSouth{grid.first_coord.dual_front(), bc_type},
                      BCNorth{grid.first_coord.dual_back(), bc_type},
                      grid.second_coord,
                      functor},
                  east_west{BCEast{grid.second_coord.dual_back(), bc_type}, BCWest{grid.second_coord.dual_front(), bc_type}, grid.first_coord, functor}
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

            void set_west_val(MatrixView &view, auto i) const
            {
                if (east_west.west.type == BoundaryCondition::BCType::first)
                {
                    view.matrix.coeffRef(view.diag) = (RealType)1.0;
                    view.matrix.coeffRef(view.neib) = (RealType)0.0;
                    view.rhs.coeffRef(0) = west_vals(i);
                    return;
                }
                else if (east_west.west.type == BoundaryCondition::BCType::second)
                {
                    view.rhs.coeffRef(0) += west_vals(i);
                    return;
                }

                assert(false && "Boundary condition is not properly set");
            }

            void set_east_val(MatrixView &view, auto i) const
            {
                if (east_west.east.type == BoundaryCondition::BCType::first)
                {
                    view.matrix.coeffRef(view.diag) = (RealType)1.0;
                    view.matrix.coeffRef(view.neib) = (RealType)0.0;
                    view.rhs.coeffRef(0ll) = east_vals(i);
                    return;
                }
                else if (east_west.east.type == BoundaryCondition::BCType::second)
                {
                    view.rhs.coeffRef(0) += east_vals(i);
                    return;
                }

                assert(false && "Boundary condition is not properly set");
            }

            void set_south_val(MatrixView &view, auto i) const
            {
                if (south_north.south.type == BoundaryCondition::BCType::first)
                {
                    view.matrix.coeffRef(view.diag) = (RealType)1.0;
                    view.matrix.coeffRef(view.neib) = (RealType)0.0;
                    view.rhs.coeffRef(0) = south_vals(i);
                    return;
                }
                else if (south_north.south.type == BoundaryCondition::BCType::second)
                {
                    view.rhs.coeffRef(0) += south_vals(i);
                    return;
                }

                assert(false && "Boundary condition is not properly set");
            }

            void set_north_val(MatrixView &view, auto i) const
            {
                if (south_north.north.type == BoundaryCondition::BCType::first)
                {
                    view.matrix.coeffRef(view.diag) = (RealType)1.0;
                    view.matrix.coeffRef(view.neib) = (RealType)0.0;
                    view.rhs.coeffRef(0ll) = north_vals(i);
                    return;
                }
                else if (south_north.north.type == BoundaryCondition::BCType::second)
                {
                    view.rhs.coeffRef(0) += north_vals(i);
                    return;
                }

                assert(false && "Boundary condition is not properly set");
            }

        protected:
            BCEastWest east_west;
            BCSouthNorth south_north;
        };
    } // BoundaryConditions
} // GPN
