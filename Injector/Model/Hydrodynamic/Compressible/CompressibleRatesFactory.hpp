#pragma once

#include <cassert>
#include <vector>

#include <Injector/Grids/Defines.h>

#include <Injector/History/InjectorRegimes.hpp>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Properties/FaceProperties.hpp>
#include <Injector/Properties/JT_FieldFactory.hpp>

#include <Injector/Solver/State2D.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <
            typename Grid2D_t, typename Well_t,
            typename History_t, typename Fluid_t,
            typename Hydrodynamics_t>
        struct CompressibleRatesFactory
        {
            using OriginalGrid = typename Grid2D_t::OriginalGrid;

            struct Solution
            {
                Solution(const auto history) 
                {
                    times.reserve(history->size());
                    states.reserve(history->size());
                }
                std::vector<RealType> times;
                std::vector<EqSolver::State::State2D> states;
            };

            Solution solution;

            CompressibleRatesFactory(
                cptr<Hydrodynamics_t> pressure_field,
                const cptr<Grid2D_t> grid2D_rocks,
                const Well_t &well,
                const ptr<History_t> history,
                const Fluid_t &fluid)
                : grid2D_rocks{grid2D_rocks},
                  grid2D{grid2D_rocks->grid2D},
                  well{well},
                  history{history},
                  fluid{fluid},
                  pressure_field{pressure_field},
                  first_size{
                      grid2D_rocks->grid2D->first_coord().mesh_size()},
                  second_size{
                      grid2D_rocks->grid2D->second_coord().mesh_size()},
                  mobility{pressure_field->mobility},
                  mobility_factor{
                      set_mobility_factor(
                          pressure_field->mobility.face_vals_axes2,
                          grid2D_rocks->first_coord().control_volumes)},
                  solution{history}
            {
            }

            /// @brief set the flow field at the next time moment
            /// @param t current time moment
            /// @param t_step step to the next time moment
            void set_flow_field(
                double t, RealType t_step)
            {
                // this method only works at FixedRate injection
                assert(history->regime() == InjectorRegimes::FixedRate);

                // the filed is updated at every time step
                // non-stationary hydrodynamics is assumed
                pressure_field->set_pressure_field(t_step, get_history_record());

                const auto mid_time{history->get_current_record().mid_time};
                if ((t < mid_time) && (t + t_step > mid_time))
                {
                    /*Properties::Pressure<Grid2D_t>*/
                    const auto &P{get_pressure_field()};
                    solution.times.push_back(t + t_step / 2.0);
                    solution.states.emplace_back(P.values());
                }

                // volumetric flow field in two directions is calculated,
                // once the pressure field is calculated
                const FaceValuesContainer axes1_value{
                    FaceValuesContainer::Zero(
                        first_size + 1ll,
                        second_size)};
                FaceValuesContainer axes2_value{
                    FaceValuesContainer::Zero(
                        first_size,
                        second_size + 1ll)};

                static_assert(Grid2D_t::l_margin == 3ll);
                // retain zero values
                //    axes2_value.leftCols(Grid2D_t::l_margin).colwise() = 0.0;
                axes2_value.col(Grid2D_t::l_margin) =
                    well.get_RFP(get_history_record());

                // face values of mobility are not defined at the domain boundaries
                assert(mobility.face_vals_axes1.rows() == first_size - 1ll);
                assert(mobility.face_vals_axes1.cols() == second_size - Grid2D_t::l_margin);
                assert(mobility.face_vals_axes2.rows() == first_size);
                const auto second_size_rock{second_size - Grid2D_t::l_margin - 1ll};
                assert(mobility.face_vals_axes2.cols() == second_size_rock);
                //    axes2_value.middleCols(Grid2D_t::l_margin + 1ll, second_size_rock) =
                //        mobility_factor;

                const auto &P{get_pressure_field()};
                for (auto col{Grid2D_t::l_margin + 1ll}, count{0ll}; count < second_size_rock; ++col, ++count)
                {
                    axes2_value.col(col) = mobility_factor.col(count) * (P.col(col-1ll) - P.col(col));
                }

                for (auto row{0ll}; row < first_size; ++row)
                {
                    for (auto col{0ll}; col < Grid2D_t::l_margin; ++col)
                        assert(axes2_value(row, col) == 0.0);
                }

                // volumetric flow field in two directions is calculated,
                // once the pressure field is calculated
                heat_flow_field =
                    std::make_shared<FaceProperties::HeatFlowField>(
                        ReservoirFlowField{
                            axes1_value,
                            axes2_value},
                        // multiply by heat capacity
                        fluid.volumetric_heat_capacity);
            }

            const auto get_spatial_JT_contribution() const
            {
                return Properties::JT_FieldFactory::create(*this).values();
            }

            const auto get_heat_flow_in_axes1() const
            {
                return heat_flow_field->axes1_as_face_normal_pos + heat_flow_field->axes1_as_face_normal_neg;
            }
            const auto get_heat_flow_in_axes2() const
            {
                return heat_flow_field->axes2_as_face_normal_pos + heat_flow_field->axes2_as_face_normal_neg;
            }

            const auto &get_heat_flow_in_axes1_pos() const
            {
                return heat_flow_field->axes1_as_face_normal_pos;
            }
            const auto &get_heat_flow_in_axes2_pos() const
            {
                return heat_flow_field->axes2_as_face_normal_pos;
            }
            const auto &get_heat_flow_in_axes1_neg() const
            {
                return heat_flow_field->axes1_as_face_normal_neg;
            }
            const auto &get_heat_flow_in_axes2_neg() const
            {
                return heat_flow_field->axes2_as_face_normal_neg;
            }
            const auto &get_pressure_field() const
            {
                return pressure_field->current_pressure();
            }

        public:
            const cptr<OriginalGrid> grid2D;
            const cptr<Grid2D_t> grid2D_rocks;
            const Well_t &well;
            const FaceProperties::Mobility<Grid2D_t> &mobility;
            const FaceValuesContainer mobility_factor;
            const ptr<History_t> history;
            const Fluid_t &fluid;

            ptr<Hydrodynamics_t> pressure_field;

        protected:
            const ptrdiff_t first_size, second_size;

            cptr<FaceProperties::HeatFlowField> heat_flow_field;
            //    cptr<FaceProperties::ReservoirFlowField> volumetric_flow_field;

            const auto get_history_record() const
            {
                return history->get_current_record();
            }

        private:
            static auto set_mobility_factor(const auto &mobility, const auto &h)
            {
                const auto pi{std::numbers::pi_v<RealType>};
                return FaceValuesContainer{mobility.colwise() * (h * 2.0 * pi)};
            }
        };

    } // FaceProperties
} // GPN
