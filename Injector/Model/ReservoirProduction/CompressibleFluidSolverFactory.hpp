#pragma once

#include <memory>
#include <vector>

#include <Injector/Properties/PhysicalField.hpp>

#include <Injector/Model/Well/WellBottomHoleRateControl.hpp>
#include <Injector/Model/Well/WellBottomHolePressureControl.hpp>
#include <Injector/Model/Well/WellBottomHoleRateControl.hpp>
#include <Injector/Model/Hydrodynamic/Compressible/CompressibleFluidSolver.hpp>

#include <Injector/Solver/FullImplicit/Solver.hpp>
#include <Injector/Solver/FullImplicit/WellRateControlSolver.hpp>

namespace GPN
{
    struct CompressibleFluidSolverFactory
    {
    protected:
        template <typename Well_t>
        const static auto make_corrected_compressibility(const Well_t &well)
        {
            using Grid2D_t = Well_t::grid_type;
            const auto &is_permeable{well.rock_field_props.base_hydrodynamics.is_permeable};
            return Properties::MediumCompressibility{
                Properties::Field<Grid2D_t>{
                    well.rock_field_props.medium_compressibility.values().colwise() + (1.0 - is_permeable.log_vals),
                    well.grid2D_rocks}};
        }
        template <typename Well_t>
        const static auto make_initial_state(const Well_t &well)
        {
            return ICFactory(well.history.start_time, well.grid2D_rocks, well.rock_field_props.base_hydrodynamics.ext_pressure);
        }

    public:
        template <
            typename History_t,
            typename Grid2D_t,
            typename CrossFlow_t,
            typename FluidField_t>
        static void create_solver(
            const Wells::BotHoleRateControl::WellBottomHoleRateControl<History_t, Grid2D_t, CrossFlow_t> &well,
            const FluidField_t &fluid_field)
        {
            using Well_t = Wells::BotHoleRateControl::WellBottomHoleRateControl<History_t, Grid2D_t, CrossFlow_t>;

            using hydro_bc_type = typename Well_t::hydro_bc_type;

            using Solver_t =
                EqSolver::FullImplicit::WellRateControlSolver<
                    Grid2D_t, Properties::MediumCompressibility<Grid2D_t>, EqSolver::EmptyConvectionField, hydro_bc_type,
                    Well_t>;

            const hydro_bc_type &bc{*(well.hydro_bc)};
            const auto initial_state{make_initial_state(well)};
            const auto ptr_rates_factory{std::make_shared<EqSolver::EmptyConvectionField>()};
            const Properties::MediumCompressibility corrected_compressibility{make_corrected_compressibility()};

            return Solver_t{
                fluid_field.face_mobility,
                well.grid2D_rocks,
                corrected_compressibility,
                ptr_rates_factory, initial_state,
                bc, well, well.history.start_time};
        }

        template <
            typename History_t,
            typename Grid2D_t,
            typename CrossFlow_t,
            typename FluidField_t>
        static void create_solver(
            const Wells::BotHolePresControl::WellBottomHolePressureControl<History_t, Grid2D_t, CrossFlow_t> &well,
            const FluidField_t &fluid_field)
        {
            using Well_t = Wells::BotHolePresControl::WellBottomHolePressureControl<History_t, Grid2D_t, CrossFlow_t>;

            using hydro_bc_type = typename Well_t::hydro_bc_type;

            using Solver_t =
                EqSolver::FullImplicit::Solver<
                    Grid2D_t, Properties::MediumCompressibility<Grid2D_t>, EqSolver::EmptyConvectionField, hydro_bc_type>;

            const hydro_bc_type &bc{*(well.hydro_bc)};
            const auto initial_state{make_initial_state(well)};
            const auto ptr_rates_factory{std::make_shared<EqSolver::EmptyConvectionField>()};
            const Properties::MediumCompressibility corrected_compressibility{make_corrected_compressibility()};

            return Solver_t{
                fluid_field.face_mobility,
                well.grid2D_rocks,
                corrected_compressibility,
                ptr_rates_factory, initial_state,
                bc, well.history.start_time};
        }
    };
} // GPN