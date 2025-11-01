#pragma once

#include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
        template <
            typename Grid2D_t, typename Fluid_t, 
            typename Well_t, typename History_t>
        struct SomeFluidField
        {
            using Axes1 = std::decay_t<typename Grid2D_t::Axes1Coordinate_t>;

            SomeFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Logs::Porosity &porosity,
                const Logs::ExternalPressure &ext_pressure,
                const Well_t &well,
                const cptr<History_t> history,
                const cptr<Grid2D_t> grid2D)
                : fluid{fluid},
                  base_hydrodynamics{base_hydrodynamics},
                  permeability{permeability},
                  porosity{porosity},
                  ext_pressure{ext_pressure},
                  P{std::make_shared<Properties::Pressure<Grid2D_t>>(
                      set_initial_pressure(ext_pressure, grid2D))},
                  P_prev{std::make_shared<Properties::Pressure<Grid2D_t>>(
                      set_initial_pressure(ext_pressure, grid2D))},
                  thickness_log{grid2D->first_coord().control_volumes},
                  well{well},
                  grid2D{grid2D},
                  history{history},
                  current_time{start_time}
            {
            }

            const Properties::Pressure<Grid2D_t> &current_pressure() const
            {
                return *P;
            }
            const auto delta_pressure() const
            {
                const auto out{(*P - *P_prev).values()};
                return out;
            }

        public:
            const Fluid_t &fluid;
            const Well_t &well;
            const cptr<Grid2D_t> grid2D;
            const ControlVolumesContainer &thickness_log;
            const Logs::Hydrodynamics::BaseHydrodynamics<Axes1>& base_hydrodynamics;
            const Logs::Permeability &permeability;
            const Logs::Porosity &porosity;
            const Logs::ExternalPressure &ext_pressure;
            const cptr<History_t> history;

        protected:
            std::shared_ptr<Properties::Pressure<Grid2D_t>> P, P_prev;
            RealType current_time;

            /// @brief Set the initial pressure replicating the pressure at the external boundary
            /// @param ext_pressure Pressure log at the external boundary
            /// @return Intial pressure field, const in r-direction
            static auto set_initial_pressure(
                const auto &ext_pressure,
                const auto grid2D)
            {
                return Properties::FieldFactory::create(ext_pressure, grid2D);
            }

            template<typename HistoryRecord_t>
            void set_time(
                const RealType time_step, 
                const HistoryRecord_t &history_record)
            {
                assert(time_step <= history_record.time_step);
                current_time += time_step;
            }
        };
    } // Hydrodynamic
} // GPN