#pragma once

#include <vector>

#include <Injector/Solver/State2D.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <
            typename Grid2D_t,
            typename Well_t,
            typename History_t,
            typename Fluid_t,
            typename Hydrodynamics_t>
        struct SomeRatesFactory
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

            SomeRatesFactory(
                cptr<Hydrodynamics_t> pressure_field,
                const cptr<Grid2D_t> grid2D_rocks,
                const ptr<Well_t> well,
                const ptr<History_t> history,
                const Fluid_t &fluid)
                : history{history},
                  pressure_field{pressure_field},
                  well{well},
                  solution{history},
                  rates{history},
                  fluid{fluid},
                  grid2D_rocks{grid2D_rocks},
                  grid2D{grid2D_rocks->grid2D},
                  first_size{
                      grid2D_rocks->grid2D->first_coord().mesh_size()},
                  second_size{
                      grid2D_rocks->grid2D->second_coord().mesh_size()}
            {
            }

            /// @brief set the flow field at the next time moment
            /// @param t current time moment
            /// @param t_step step to the next time moment
            template <typename Functor_t>
            void set_flow_field(
                double t, RealType t_step, Functor_t &&collector_rates)
            {
                // this method only works at FixedRate injection
                assert((history->regime() == InjectorRegimes::FixedRate) ||
                       (history->regime() == InjectorRegimes::FixedBottomHolePressure));

                // the filed is updated at every time step
                // non-stationary hydrodynamics is assumed
                pressure_field->set_pressure_field(t_step, get_history_record());
                // set the flow in all cells of the well,
                // taking RFP and WFP into account
                well->set_well_flow_field(get_history_record(), pressure_field->get_rock_pressure());

                const auto &P{get_pressure_field()};
                
                for (auto col{Grid2D_t::l_margin + 1ll}; col < second_size; ++col)
                {
                    for (auto row{0ll}; row < P.rows(); ++row)
                    {
                        assert(!std::isnan(P.value(row, col - 1ll)) && !std::isinf(P.value(row, col - 1ll)));
                        assert(!std::isnan(P.value(row, col)) && !std::isinf(P.value(row, col)));
                    }
                }

                const auto mid_time{history->get_current_record().mid_time};
                if ((t < mid_time) && (t + t_step > mid_time))
                {
                    /*Properties::Pressure<Grid2D_t>*/
                    solution.times.push_back(t + t_step / 2.0);
                    solution.states.emplace_back(P.values());
                }

                // volumetric heat flow field in two directions is calculated,
                // once the pressure field is calculated
                FaceValuesContainer axes1_value{
                    FaceValuesContainer::Zero(
                        first_size + 1ll,
                        second_size)};
                axes1_value.leftCols(Grid2D_t::l_margin) = well->flow_axes1_value;
                // other verticle flow components are zero

                for (auto row{0ll}; row < axes1_value.leftCols(Grid2D_t::l_margin).rows(); ++row)
                    for (auto col{0ll}; col < axes1_value.leftCols(Grid2D_t::l_margin).cols(); ++col)
                        assert(!std::isnan(axes1_value(row, col)) && !std::isinf(axes1_value(row, col)));

                FaceValuesContainer axes2_value{
                    FaceValuesContainer::Zero(
                        first_size,
                        second_size + 1ll)};
                axes2_value.leftCols(Grid2D_t::l_margin + 1ll) = well->flow_axes2_value;

                for (auto row{0ll}; row < axes2_value.leftCols(Grid2D_t::l_margin + 1ll).rows(); ++row)
                    for (auto col{0ll}; col < axes2_value.leftCols(Grid2D_t::l_margin + 1ll).cols(); ++col)
                        assert(!std::isnan(axes2_value(row, col)) && !std::isinf(axes2_value(row, col)));

                for (auto col{Grid2D_t::l_margin + 1ll}, count{0ll}; col < second_size; ++col, ++count)
                    axes2_value.col(col) = collector_rates(count, col, P); // = well->rfp

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

            const auto &get_pressure_field() const
            {
                return pressure_field->current_pressure();
            }

            const auto get_spatial_JT_contribution() const
            {
                return Properties::JT_FieldFactory::create_spatial(*this).values();
            }

            const auto get_temporal_JT_contribution() const
            {
                return Properties::JT_FieldFactory::create_temporal(*this).values();
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
                const auto &values{heat_flow_field->axes1_as_face_normal_pos};
                for (auto row{0ll}; row < values.rows(); ++row)
                    for (auto col{0ll}; col < values.cols(); ++col)
                        assert(!std::isnan(values(row, col)) && !std::isinf(values(row, col)));
                return values;
            }
            const auto &get_heat_flow_in_axes2_pos() const
            {
                const auto &values{heat_flow_field->axes2_as_face_normal_pos};
                for (auto row{0ll}; row < values.rows(); ++row)
                    for (auto col{0ll}; col < values.cols(); ++col)
                        assert(!std::isnan(values(row, col)) && !std::isinf(values(row, col)));
                return values;
            }
            const auto &get_heat_flow_in_axes1_neg() const
            {
                const auto &values{heat_flow_field->axes1_as_face_normal_neg};
                for (auto row{0ll}; row < values.rows(); ++row)
                    for (auto col{0ll}; col < values.cols(); ++col)
                        assert(!std::isnan(values(row, col)) && !std::isinf(values(row, col)));
                return values;
            }
            const auto &get_heat_flow_in_axes2_neg() const
            {
                const auto &values{heat_flow_field->axes2_as_face_normal_neg};
                for (auto row{0ll}; row < values.rows(); ++row)
                    for (auto col{0ll}; col < values.cols(); ++col)
                        assert(!std::isnan(values(row, col)) && !std::isinf(values(row, col)));
                return values;
            }

        public:
            const ptr<History_t> history;
            const Fluid_t &fluid;
            ptr<Hydrodynamics_t> pressure_field;
            const ptr<Well_t> well;
            const cptr<Grid2D_t> grid2D_rocks;
            const cptr<OriginalGrid> grid2D;

            Solution solution;
            Solution rates;

        protected:
            const ptrdiff_t first_size, second_size;

            cptr<FaceProperties::HeatFlowField> heat_flow_field;
            cptr<FaceProperties::ReservoirFlowField> volumetric_flow_field;

            const auto get_history_record() const
            {
                return history->get_current_record();
            }
        };
    } // Hydrodynamic
} // GPN