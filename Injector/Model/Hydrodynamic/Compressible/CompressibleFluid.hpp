#pragma once

// #include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/Model/Hydrodynamic/SomeFluidField.hpp>


namespace GPN
{
    namespace Hydrodynamic
    {
        template<typename Grid2D_t, typename Fluid_t, typename Well_t>
        struct CompressibleFluid
            : public SomeFluidField<Grid2D_t, Fluid_t, Well_t>
        {
            using SomeFluidField<Grid2D_t, Fluid_t, Well_t>::ext_pressure;
            using SomeFluidField<Grid2D_t, Fluid_t, Well_t>::P;
            using SomeFluidField<Grid2D_t, Fluid_t, Well_t>::grid2D;

            CompressibleFluid(
                const RealType start_time,
                const Fluid_t& fluid,
                const Logs::Permeability &permeability,
                const Logs::ExternalPressure &ext_pressure,
                const Well_t &well,
                const cptr<Grid2D_t> grid2D)
                : SomeFluidField<Grid2D_t, Fluid_t, Well_t>{
                      start_time, fluid,
                      permeability,
                      ext_pressure,
                      well, grid2D}
            {
            }

            template<typename HistoryRecord_t>
            void set_pressure_field(
                const RealType time_step, 
                const HistoryRecord_t &history_record)
            {
                // calculate current pressure well and cement-sandwich
                auto v{CellNodesContainer2D::Zero(
                    grid2D->first_coord.mesh_size(), 
                    grid2D->second_coord.mesh_size())};
                v.colwise() = ext_pressure.log_vals;

                P = std::make_shared<Properties::Pressure<Grid2D_t>>(
                    v,
                    grid2D);
            }            
        };
    } // Hydrodynamic
} // GPN