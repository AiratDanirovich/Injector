#pragma once

#include <numbers>
#include <cassert>

#include <Injector/Model/Hydrodynamic/SomeFluidField.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
        template <typename Grid2D_t, typename Fluid_t, typename Well_t>
        struct IncompressibleFluidField
            : public SomeFluidField<Grid2D_t, Fluid_t, Well_t>
        {
            using SomeFluidField<Grid2D_t, Fluid_t, Well_t>::ext_pressure;
            using SomeFluidField<Grid2D_t, Fluid_t, Well_t>::P;
            using SomeFluidField<Grid2D_t, Fluid_t, Well_t>::grid2D;
            using SomeFluidField<Grid2D_t, Fluid_t, Well_t>::well;

            IncompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Logs::ExternalPressure &ext_pressure,
                const Well_t &well,
                const cptr<Grid2D_t> grid2D)
                : SomeFluidField<Grid2D_t, Fluid_t, Well_t>{
                      start_time, fluid,
                      permeability,
                      ext_pressure,
                      well, grid2D},
                  auxillary_term{set_auxillary_term(fluid, ext_pressure, permeability, grid2D)}
            {
            }

            template<typename HistoryRecord_t>
            void set_pressure_field(
                const RealType time_step, 
                const HistoryRecord_t &history_record)
            {
                const StepPropertyContainer RFP{well.get_RFP(history_record)};
                // calculate current pressure well and cement-sandwich
                auto v{((auxillary_term.colwise() * RFP).colwise() + ext_pressure.log_vals).eval()};
                P = std::make_shared<Properties::Pressure<Grid2D_t>>(
                    v,
                    grid2D);

                SomeFluidField<Grid2D_t, Fluid_t, Well_t>::set_time(time_step, history_record);
            }

            const CellNodesContainer2D auxillary_term;

        private:
            static CellNodesContainer2D set_auxillary_term(
                const auto &fluid, const auto &ext_pressure, const auto &permeability, const auto grid2D)
            {
                const auto thickness_log{grid2D->first_coord.control_volumes};
                const auto &r_grid{grid2D->second_coord};
                const auto r_max{r_grid.dual_back()};
                const auto pi{std::numbers::pi_v<RealType>};

                auto r_nodes{r_grid.mesh_nodes};
                // account for the well,
                // the pressure is contant within the sandface radius
                r_nodes.head(3ll) = r_grid.dual_nodes(3ll);

                const StepPropertyContainer temp1{(-fluid.viscosity) * (r_nodes / r_max).log()};
                StepPropertyContainer temp2{1.0 / ((2.0 * pi) * thickness_log * permeability.log_vals)};

                // account for layers with zero permeability
                for (auto row{0ll}; row < permeability.size(); ++row)
                    if (permeability(row) == 0.0)
                        temp2(row) = 0.0;

                const CellNodesContainer2D temp{(temp2.matrix() * temp1.transpose().matrix()).array()};

                assert(temp.rows() == grid2D->first_coord.mesh_size());
                assert(temp.cols() == grid2D->second_coord.mesh_size());

                return temp;
            }
        };
    } // Hydrodynamic
} // GPN