#pragma once

#include <numbers>

#include <Injector/Grids/Defines.h>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Hydrodynamic/SomeFluidField.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluidSolver.hpp>

#include <Injector/Solver/FullImplicit/Solver.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
        template <
            typename Grid2D_t, typename Fluid_t,
            typename Well_t, typename History_t>
        struct CompressibleFluidField
            : public SomeFluidField<typename Grid2D_t::OriginalGrid, Fluid_t, Well_t, History_t>
        {
            using OriginalGrid = typename Grid2D_t::OriginalGrid;
            using Base = SomeFluidField<OriginalGrid, Fluid_t, Well_t, History_t>;
            using Base::ext_pressure;
            using Base::grid2D;
            using Base::P;
            using Base::well;

            using Solver_t =
                EqSolver::FullImplicit::Solver<
                    Grid2D_t,
                    Properties::MediumCompressibility<Grid2D_t>,
                    EqSolver::EmptyConvectionField,
                    HydroBC>;

            CompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                const Logs::ExternalPressure &ext_pressure,
                const Well_t &well,
                const cptr<History_t> history,
                const cptr<Grid2D_t> grid2D_rocks)
                : Base{
                      start_time, fluid,
                      permeability,
                      ext_pressure,
                      well, history, grid2D_rocks->grid2D},
                  mobility{FaceProperties::Rocks::RocksFaceProps{
                      rock_field_props,
                      grid2D_rocks}.mobility}, 
                  first_size{
                    grid2D_rocks->grid2D->first_coord().mesh_size()}, 
                  second_size{
                    grid2D_rocks->grid2D->second_coord().mesh_size()},
                  inv_mobility{set_inv_mobility(rock_field_props, well, grid2D_rocks)}, 
                  grid2D_rocks{grid2D_rocks}
            {
                solver =
                    std::make_unique<CompressibleFluidSolver<Solver_t>>(set_solver(
                        start_time, history, 
                        rock_field_props, 
                        ext_pressure, 
                        well.RFP_weights,
                        mobility,
                        grid2D_rocks));
            }

            template <typename HistoryRecord_t>
            void set_pressure_field(
                const RealType time_step,
                const HistoryRecord_t &record)
            {
                solver->advance(time_step);
                const GridNodeValues2D &rock_P{solver->get_state().cur_state};
                GridNodeValues2D out{GridNodeValues2D::Zero(first_size, second_size)};
                out.rightCols(second_size - Grid2D_t::l_margin) = rock_P;
                out.leftCols(Grid2D_t::l_margin).colwise() = rock_P.col(0ll) + record.rate * inv_mobility;

                P = std::make_shared<Properties::Pressure<OriginalGrid>>(
                    std::move(out),
                    grid2D);
            }

            const auto& get_mobility() const
            {
                return solver->mobility;
            }

            const FaceProperties::Mobility<Grid2D_t> mobility;
        private:
            const Eigen::ArrayX<RealType> inv_mobility;
            std::unique_ptr<CompressibleFluidSolver<Solver_t>> solver;
            const ptrdiff_t first_size, second_size;
            const cptr<Grid2D_t> grid2D_rocks;

            static auto set_inv_mobility(
                const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                const Well_t &well,
                cptr<Grid2D_t> grid2D_rocks)
            {
                const auto r3{grid2D_rocks->second_coord().mesh_nodes(0ll)};
                const auto r2_face{grid2D_rocks->second_coord().dual_nodes(0ll)};
                const auto two_pi{2.0 * std::numbers::pi_v<RealType>};
                const auto is_permeable{rock_field_props.base_hydrodynamics.is_permeable.log_vals};

                Eigen::ArrayX<RealType> out{well.RFP_weights / (
                    two_pi / std::log(r3 / r2_face) * 
                    rock_field_props.mobility_axes2.col(Grid2D_t::l_margin) * 
                    grid2D_rocks->first_coord().volumes())};

                    for(auto row{0ll}; row < grid2D_rocks->first_coord().mesh_size(); ++row)
                    {
                        if(is_permeable(row) == 0.0)
                            out(row) = 0.0;
                    }
                return out;
            }

            static auto set_solver(
                const RealType start_time,
                const cptr<History_t> history,
                const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                const Logs::ExternalPressure &ext_pressure,
                const auto rfp,
                const FaceProperties::Mobility<Grid2D_t>& mobility,
                const cptr<Grid2D_t> grid2D_rocks)
            {
                // FaceProperties::Rocks::RocksFaceProps
                //     rock_face_props{
                //         rock_field_props,
                //         grid2D_rocks};

                const HydroBC bc{
                    grid2D_rocks,
                    std::make_shared<const FunctorBC<
                        Well_t, History_t, Grid2D_t>>(
                        history, ext_pressure, rfp, grid2D_rocks)};

                const auto initial_state{ICFactory(start_time, grid2D_rocks, ext_pressure)};

                const auto ptr_rates_factory{std::make_shared<EqSolver::EmptyConvectionField>()};

                const auto &is_permeable{rock_field_props.base_hydrodynamics.is_permeable};

                const Properties::MediumCompressibility corrected_compressibility{
                    Properties::Field<Grid2D_t>{
                        rock_field_props.medium_compressibility.values().colwise() + (1.0 - is_permeable.log_vals),
                        grid2D_rocks}};

                using Solver_t = decltype(EqSolver::FullImplicit::Solver{
                    /*rock_face_props.*/mobility,
                    grid2D_rocks,
                    corrected_compressibility,
                    ptr_rates_factory, initial_state,
                    bc, start_time});

                auto solver_ptr{std::make_shared<Solver_t>(
                    /*rock_face_props.*/mobility,
                    grid2D_rocks,
                    corrected_compressibility,
                    ptr_rates_factory, initial_state,
                    bc, start_time)};

                return solver_ptr;
            }
        };
    } // Hydrodynamic
} // GPN