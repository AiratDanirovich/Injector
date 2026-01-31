#pragma once

#include <numbers>
#include <memory>

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/Logs.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Hydrodynamic/SomeFluidField.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluidSolver.hpp>
#include <Injector/Model/ReservoirProduction/CompressibleFluidSolverFactory.hpp>

#include <Injector/Solver/State2D.hpp>
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
            using Base::P_prev;
            using Base::well;
            using hydro_bc_type = typename Well_t::hydro_bc_type;

            using Solver_t =
                typename Well_t:: template solver_type<
                    Grid2D_t,
                    Properties::MediumCompressibility<Grid2D_t>,
                    EqSolver::EmptyConvectionField,
                    hydro_bc_type>;

            CompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Properties::Rocks::RocksProps<Grid2D_t> &
                    rock_field_props,
                const Well_t &well,
                const cptr<History_t> history,
                const cptr<Grid2D_t> grid2D_rocks)
                : Base{
                      start_time, fluid,
                      rock_field_props.base_hydrodynamics.permeability,
                      rock_field_props.base_hydrodynamics.porosity,
                      rock_field_props.base_hydrodynamics.ext_pressure,
                      well, history, grid2D_rocks->grid2D},
                  face_mobility{
                    FaceProperties::Rocks::RocksFaceProps{rock_field_props, grid2D_rocks}.mobility},
                  first_size{grid2D_rocks->grid2D->first_coord().mesh_size()},
                  second_size{grid2D_rocks->grid2D->second_coord().mesh_size()},
                  grid2D_rocks{grid2D_rocks}
            {
                solver = std::make_unique<CompressibleFluidSolver<Solver_t>>(
                        CompressibleFluidSolverFactory::create_solver_ptr(well, *this)
                );
                // solver =
                //     std::make_unique<CompressibleFluidSolver<Solver_t>>(set_solver(
                //         start_time, history,
                //         rock_field_props,
                //         mobility,
                //         grid2D_rocks));
            }

            template <typename HistoryRecord_t>
            void set_pressure_field(
                const RealType time_step,
                const HistoryRecord_t &history_record)
            {
                solver->advance(time_step);
                const GridNodeValues2D &rock_P{get_rock_pressure()};
                GridNodeValues2D out{GridNodeValues2D::Zero(first_size, second_size)};
                // set pressure in reservoir as a solution of respective problem
                out.rightCols(second_size - Grid2D_t::l_margin) = rock_P;
                // get pressure at the sandface
                out.leftCols(Grid2D_t::l_margin).colwise() = 
                    well.get_pressure_at_sandface(
                        history_record, rock_P);

                for (auto row{0ll}; row < out.rows(); ++row)
                    for (auto col{0ll}; col < out.cols(); ++col)
                        assert(!std::isnan(out(row, col)) && !std::isinf(out(row, col)));

                P_prev = P;

                P = std::make_shared<Properties::Pressure<OriginalGrid>>(
                    std::move(out),
                    grid2D);

                Base::set_time(time_step, history_record);
            }

            const auto &get_mobility() const
            {
                return solver->mobility;
            }
            const auto &get_rock_pressure() const
            {
                return solver->get_state().cur_state;
            }

            const auto& get_solver() const
            {
                return solver;
            }

            const FaceProperties::Mobility<Grid2D_t> face_mobility;

        private:
            std::unique_ptr<CompressibleFluidSolver<Solver_t>> solver;
            const ptrdiff_t first_size, second_size;
            const cptr<Grid2D_t> grid2D_rocks;
        };
    } // Hydrodynamic
} // GPN