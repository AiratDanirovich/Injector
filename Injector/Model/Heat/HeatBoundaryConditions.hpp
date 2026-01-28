#pragma once

#include <memory>

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/History/History.hpp>

namespace GPN
{
    namespace Heat
    {
        /// @brief Boundary condition for temperature equation. 
        /// Depends on the flow hydrodynamic model Fieldfactory_t in reservoir.
        /// @tparam Fieldfactory_t FaceProperties::CompressibleRatesFactory or FaceProperties::IncompressibleRatesFactory
        template <typename Fieldfactory_t>
        struct HeatBC : public BoundaryConditions::GeneralBC
        {
            template <typename Grid2D_t>
            HeatBC(const cptr<Grid2D_t> &grid,
                   cptr<const BCFunctorBase> functor,
                   const cptr<Fieldfactory_t> field_factory)
                : BoundaryConditions::GeneralBC{grid, functor, BoundaryCondition::BCType::second},
                field_factory{field_factory},
                history{field_factory->history}
            {
            }

            void set_bc_type(const RealType t)
            {
                this->t = t;
            //    bc_types[south_id] = BoundaryCondition::BCType::second; // top boundary
            //    bc_types[north_id] = BoundaryCondition::BCType::second; // bottom boundary
            //    bc_types[west_id] = BoundaryCondition::BCType::second;  // well axis of symmetry
            if(history->rate() > 0.0)
                bc_types[east_id] = BoundaryCondition::BCType::second;  // external contour
            else
                bc_types[east_id] = BoundaryCondition::BCType::first;  // external contour
            }

        protected:
            const cptr<Fieldfactory_t> field_factory;
            const cptr<History> history;

            bool is_outward_flux() const
            {
                
            }
        };
    } // Heat
} // GPN