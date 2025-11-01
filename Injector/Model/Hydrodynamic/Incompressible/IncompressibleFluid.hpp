#pragma once

#include <numbers>
#include <cassert>

#include <Injector/Properties/Logs.hpp>
#include <Injector/Model/Hydrodynamic/SomeFluidField.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
        template <typename Grid2D_t, typename Fluid_t, typename Well_t, typename History_t>
        struct IncompressibleFluidField
            : public SomeFluidField<typename Grid2D_t::OriginalGrid, Fluid_t, Well_t, History_t>
        {
            using OriginalGrid = typename Grid2D_t::OriginalGrid;
            using Base = SomeFluidField<OriginalGrid, Fluid_t, Well_t, History_t>;
            using Base::ext_pressure;
            using Base::P;
            using Base::grid2D;
            using Base::well;

            IncompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability, // delete
                const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                const Logs::ExternalPressure &ext_pressure, // delete
                const Well_t &well,
                const cptr<History_t> history,
                const cptr<Grid2D_t> grid2D_rocks)
                : Base{
                      start_time, fluid,
                      rock_field_props.base_hydrodynamics.permeability,
                      rock_field_props.base_hydrodynamics.ext_pressure,
                      well, history, grid2D_rocks->grid2D},
                  auxillary_term{set_auxillary_term(
                    fluid, rock_field_props.base_hydrodynamics.ext_pressure, 
                    rock_field_props.base_hydrodynamics.permeability, 
                    grid2D_rocks)},
                  first_size{
                    grid2D_rocks->grid2D->first_coord().mesh_size()}, 
                  second_size{
                    grid2D_rocks->grid2D->second_coord().mesh_size()}
            {
            }

            template<typename HistoryRecord_t>
            void set_pressure_field(
                const RealType time_step, 
                const HistoryRecord_t &history_record)
            {
                const StepPropertyContainer RFP{well.get_RFP(history_record)};
                // calculate current pressure in well and cement-sandwich
                const auto rock_P{((auxillary_term.colwise() * RFP).colwise() + ext_pressure.log_vals).eval()};

                GridNodeValues2D out{GridNodeValues2D::Zero(first_size, second_size)};
                out.rightCols(second_size - Grid2D_t::l_margin+1ll) = rock_P;

                const auto left_col{rock_P.col(0ll)};
                out.leftCols(Grid2D_t::l_margin).colwise() = rock_P.col(0ll);

                P = std::make_shared<Properties::Pressure<OriginalGrid>>(
                    std::move(out),
                    grid2D);

                Base::set_time(time_step, history_record);
            }

        private:
            const CellNodesContainer2D auxillary_term;
            const ptrdiff_t first_size, second_size;

            static CellNodesContainer2D set_auxillary_term(
                const auto &fluid, const auto &ext_pressure, 
                const auto &permeability, const auto grid2D_rocks)
            {
                const auto thickness{grid2D_rocks->first_coord().volumes()};
                const auto &r_grid{grid2D_rocks->second_coord()};
                const auto r_max{r_grid.dual_back()};
                const auto r_min{r_grid.dual_front()}; // sandface
                const auto pi{std::numbers::pi_v<RealType>};

                MeshNodesContainerT r_nodes{MeshNodesContainerT::Zero(r_grid.dual_size())};
                r_nodes(0ll) = r_min;
                r_nodes.rightCols(r_grid.mesh_size()) = r_grid.mesh_nodes;

                const StepPropertyContainer temp1{(-fluid.viscosity) * (r_nodes / r_max).log()};
                StepPropertyContainer temp2{1.0 / ((2.0 * pi) * thickness * permeability.log_vals)};

                // account for layers with zero permeability
                for (auto row{0ll}; row < permeability.size(); ++row)
                    if (permeability(row) == 0.0)
                        temp2(row) = 0.0;

                const CellNodesContainer2D temp{(temp2.matrix() * temp1.transpose().matrix()).array()};

                assert(temp.rows() == grid2D_rocks->first_coord().mesh_size());
                assert(temp.cols() == grid2D_rocks->second_coord().mesh_size());

                return temp;
            }
        };
    } // Hydrodynamic
} // GPN