#pragma once

#include <memory>
#include <vector>
#include <array>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids2D.hpp>

namespace GPN
{
    namespace BoundaryConditions
    {
        struct GeneralBC
        {
            struct BoundaryCondition
            {
                // sides of computational domain
                enum Side
                {
                    East, West,
                    North, South
                };

                // types of boundary conditions
                enum BCType
                {
                    undef,
                    first,
                    second,
                    third
                };
            };
            
            struct BCFunctorBase
            {
                using BCType = BoundaryConditions::GeneralBC::BoundaryCondition::BCType;
                virtual RealType operator()(
                    const ptrdiff_t x, RealType y, const RealType t,
                    const BoundaryCondition::BCType bc_type = BoundaryCondition::BCType::second) const = 0;
                virtual RealType operator()(
                    RealType x, const ptrdiff_t y, const RealType t,
                    const BoundaryCondition::BCType bc_type = BoundaryCondition::BCType::second) const = 0;
            };

            template <typename Grid2D_t>
            GeneralBC(
                const cptr<Grid2D_t> &grid,
                const cptr<const BCFunctorBase> functor,
                BoundaryCondition::BCType bc_type = BoundaryCondition::undef)
                : 
            GeneralBC(
                grid, functor,
                std::array<BoundaryCondition::BCType, 4ull>{bc_type,bc_type,bc_type,bc_type})
            {
            }
            
            template <typename Grid2D_t>
            GeneralBC(
                const cptr<Grid2D_t> &grid,
                const cptr<const BCFunctorBase> functor,
                std::array<BoundaryCondition::BCType, 4ull> bc_types
            )
                :   bc_types{bc_types},
                    fixed_coords{
                        grid->first_coord().dual_front(),
                        grid->first_coord().dual_back(),
                        grid->second_coord().dual_front(), 
                        grid->second_coord().dual_back()},
                    functor{functor},
                    t{0.0}
            {
            }

            GeneralBC(const GeneralBC &) = default;

            void set_bc_type(const RealType t) = delete;

            template<typename MatrixView_t>
            void set_west_val(MatrixView_t &view, const auto i) const
            {
                const RealType y{fixed_coords[west_id]};
                const BoundaryCondition::BCType 
                    bc_type{bc_types[west_id]};
                if (bc_type == BoundaryCondition::BCType::first)
                {
                    view.set_type_I((*functor)(i, y, t, bc_type));
                    return;
                }
                else if (bc_type == BoundaryCondition::BCType::second)
                {
                    view.add_rhs_type_II((*functor)(i, y, t, bc_type));
                    return;
                }

                assert(false && "West boundary condition is not properly set");
            }

            template<typename MatrixView_t>
            void set_east_val(MatrixView_t &view, const auto i) const
            {
                const RealType y{fixed_coords[east_id]};
                const BoundaryCondition::BCType 
                    bc_type{bc_types[east_id]};
                if (bc_type == BoundaryCondition::BCType::first)
                {
                    view.set_type_I((*functor)(i, y, t, bc_type));
                    return;
                }
                else if (bc_type == BoundaryCondition::BCType::second)
                {
                    view.add_rhs_type_II((*functor)(i, y, t, bc_type));
                    return;
                }

                assert(false && "East boundary condition is not properly set");
            }

            template<typename MatrixView_t>
            void set_south_val(MatrixView_t &view, const auto i) const
            {
                const RealType x{fixed_coords[south_id]};
                const BoundaryCondition::BCType 
                    bc_type{bc_types[south_id]};
                if (bc_type == BoundaryCondition::BCType::first)
                {
                    view.set_type_I((*functor)(x, i, t, bc_type));
                    return;
                }
                else if (bc_type == BoundaryCondition::BCType::second)
                {
                    view.add_rhs_type_II((*functor)(x, i, t, bc_type));
                    return;
                }

                assert(false && "South boundary condition is not properly set");
            }

            template<typename MatrixView_t>
            void set_north_val(MatrixView_t &view, const auto i) const
            {
                const RealType x{fixed_coords[north_id]};
                const BoundaryCondition::BCType 
                    bc_type{bc_types[north_id]};
                if (bc_type == BoundaryCondition::BCType::first)
                {
                    view.set_type_I((*functor)(x, i, t, bc_type));
                    return;
                }
                else if (bc_type == BoundaryCondition::BCType::second)
                {
                    view.add_rhs_type_II((*functor)(x, i, t, bc_type));
                    return;
                }

                assert(false && "North boundary condition is not properly set");
            }

        protected:
            RealType t;
            std::array<BoundaryCondition::BCType, 4ull> bc_types;
            const std::array<RealType, 4ull> fixed_coords;
            const cptr<const BCFunctorBase> functor;

            const size_t west_id{2ull}, east_id{3ull}, south_id{0ull}, north_id{1ull};
        };
    } // BoundaryConditions
} // GPN
