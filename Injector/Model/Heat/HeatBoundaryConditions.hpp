#pragma once

#include <memory>

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/History/History.hpp>

namespace GPN
{
    namespace Heat
    {
        template <typename Fieldfactory_t
        /*FaceProperties::CompressibleRatesFactory, 
        FaceProperties::IncompressibleRatesFactory*/
        >
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
                // external contour
                const auto ext_collector_pressure{field_factory->get_pressure_field().values().rightCols(2ll)};
                for (auto row{0ll}; row < field_factory->grid2D->first_coord().mesh_size(); ++row)
                {
                    if (ext_collector_pressure(row, 1) < ext_collector_pressure(row,0))
                        // liquid flows OUT of the well
                        east_bc_type[row] = BoundaryCondition::BCType::second; 
                    else
                        // liquid flows from infinity TOWARDS the well
                        east_bc_type[row] = BoundaryCondition::BCType::first;
                }
            }

        protected:
            const cptr<Fieldfactory_t> field_factory;
            const cptr<History> history;
        };
    } // Heat
} // GPN