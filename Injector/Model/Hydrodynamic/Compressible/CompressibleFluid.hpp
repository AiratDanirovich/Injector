#pragma once

#include <Injector/Grids/Defines.h>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Hydrodynamic/SomeFluidField.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluidSolver.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
        template <
        typename Grid2D_t, typename Fluid_t, 
        typename Well_t, typename History_t>
        struct CompressibleFluidField
            : public SomeFluidField<Grid2D_t, Fluid_t, Well_t, History_t>
        {
            using Base = SomeFluidField<Grid2D_t, Fluid_t, Well_t, History_t>;
            using Base::ext_pressure;
            using Base::P;
            using Base::grid2D;
            using Base::well;

            CompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                const Logs::ExternalPressure &ext_pressure,
                const Well_t &well,
                const cptr<History_t> history,
                const cptr<Grid2D_t> grid2D)
                : SomeFluidField<Grid2D_t, Fluid_t, Well_t, History_t>{
                      start_time, fluid,
                      permeability,
                      ext_pressure,
                      well, history, grid2D}
            {
                // const HydroBC bc{
                //     grid2D,
                //     std::make_shared<const FunctorBC<
                //         Well_CrossFlow,
                //         History>>(
                //         history, ext_pressure, well.RFP_weights, grid2D)};
            //    CompressibleFluidSolver solver{ext_pressure, grid2D};
            }

            template <typename HistoryRecord_t>
            void set_pressure_field(
                const RealType time_step,
                const HistoryRecord_t &history_record)
            {
                // calculate current pressure well and cement-sandwich
                auto v{
                    ext_pressure.log_vals.replicate(
                        1ll, grid2D->second_coord().mesh_size())};

                P = std::make_shared<Properties::Pressure<Grid2D_t>>(
                    v,
                    grid2D);
            }
        };
    } // Hydrodynamic
} // GPN