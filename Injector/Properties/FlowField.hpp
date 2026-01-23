#pragma once
#include <numeric>
#include <algorithm>
#include <cmath>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/LogsFactory.hpp>

namespace GPN
{
    namespace Logs
    {
        /// @brief Generates the z-component of the flow field as a function of r (i.e., other coordinate)
        struct ZFlowRateLog
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            ZFlowRateLog(
                const StepPropertyGrid &data)
                : StepPropertyGrid{data},
                  AssertNonNegative{data}
            {
            }
        };

        struct ZFlowRateLogFactory
        {
            template <typename Grid_t>
            static auto create(
                RealType well_rate,
                const Grid_t &grid1D)
            {
                assert(!std::isnan(well_rate) && !std::isinf(well_rate));
                const auto temp{make_rates(well_rate, grid1D)};
                return ZFlowRateLog{temp};
            }

        private:
            template <typename Grid_t>
            static StepPropertyGrid make_rates(
                RealType well_rate,
                const Grid_t &grid1D)
            {
                assert(!std::isnan(well_rate) && !std::isinf(well_rate));

                StepPropertyContainer verticle_rates_vals{
                    StepPropertyContainer::Zero(grid1D.mesh_size())};
                verticle_rates_vals(0ll) = well_rate;

                return {verticle_rates_vals,
                        grid1D};
            }
        };
    } // Logs

    namespace FaceProperties
    {
        struct FlowFieldFactory
        {
            template <typename Grid2D_t>
            static auto flow_in_dir1(
                const Logs::StepPropertyGrid &log,
                const Grid2D_t &grid2D)
            {
                assert(grid2D.second_coord().mesh_size() == log.size());
                FaceValuesContainer face_vals(grid2D.first_coord().dual_size(), log.size());
                face_vals.rowwise() = log.log_vals.transpose();
                return face_vals;
            }

            template <typename Grid2D_t>
            static auto flow_in_dir2(
                const StepPropertyContainer &log,
                const Grid2D_t &grid2D)
            {
                assert(grid2D.first_coord().mesh_size() == log.size());
                FaceValuesContainer face_vals(log.size(), grid2D.second_coord().dual_size());
                face_vals.colwise() = log;
                return face_vals;
            }

            template <typename Grid2D_t>
            static auto flow_in_dir2(
                const Logs::StepPropertyGrid &log,
                const Grid2D_t &grid2D)
            {
                return flow_in_dir2(log.log_vals, grid2D);
            }
        };

        struct ReservoirFlowField
        {
            ReservoirFlowField(ReservoirFlowField &&) noexcept = default;

            ReservoirFlowField(
                const FaceValuesContainer &axes1_value,
                const FaceValuesContainer &axes2_value)
                : axes1_as_face_normal{axes1_value},
                  axes2_as_face_normal{axes2_value}
            {
            }

            template <typename Grid2D_t>
            ReservoirFlowField(
                const Logs::StepPropertyGrid &axes1_value,
                const Logs::StepPropertyGrid &axes2_value,
                const Grid2D_t &grid2D)
                : ReservoirFlowField{
                      FlowFieldFactory::flow_in_dir1(axes1_value, grid2D),
                      FlowFieldFactory::flow_in_dir2(axes2_value, grid2D)}
            {
            }

            FaceValuesContainer axes1_as_face_normal;
            FaceValuesContainer axes2_as_face_normal;
        };

        void multiply(ReservoirFlowField &flow, RealType fluid_vol_heatcapacity)
        {
            flow.axes1_as_face_normal *= fluid_vol_heatcapacity;
            flow.axes2_as_face_normal *= fluid_vol_heatcapacity;
        }

        struct FlowFactory
        {
            template <typename Record_t, typename Well_t, typename Grid2D_t>
            static auto create_from_well(
                const Record_t &history_record,
                const Well_t &well,
                const Grid2D_t &grid2D)
            {
                // flow at the tube radius
                const auto wfp_vals{well.get_WFP(history_record)};
                const auto rfp_vals{well.get_RFP(history_record)};
                const auto cement_verticle_flow_vals{well.get_verticle_cement_flow(history_record)};
#pragma region AXES2-AS-FACENORMAL
                // assume reservoir flow profile by default in the entire domain
                auto axes2_as_face_normal{FlowFieldFactory::flow_in_dir2(rfp_vals, grid2D)};
                // fix horizontal flow to take boundary conditions, the well and its competion into account
                axes2_as_face_normal.col(0ll) = 0.0;      // boundary condition, zero flux at the axis of symmetry
                axes2_as_face_normal.col(1ll) = wfp_vals; // flow at the tube inner radius
                axes2_as_face_normal.col(2ll) = wfp_vals; // flow at the cement inner radius
#pragma endregion
#pragma region AXES1-AS-FACENORMAL
                // assume by default zero verticle flux in the entire domain
                FaceValuesContainer axes1_as_face_normal{
                    FaceValuesContainer::Zero(
                        grid2D.first_coord().dual_size(),
                        grid2D.second_coord().mesh_size())};

                { // set the flow in tube
                    // cumsum of rfp flow rates
                    LogValuesContainer cum_sum{LogValuesContainer::Zero(grid2D.first_coord().dual_size())};
                    std::partial_sum(wfp_vals.cbegin(), wfp_vals.cend(), cum_sum.begin() + 1ll, std::plus<RealType>());
                    // residual flowrate along the well
                    const LogValuesContainer z_flow{history_record.rate - cum_sum};

                    axes1_as_face_normal.col(0ll) = z_flow;
                }
                {
                    // the verticle flow in tube-wall::annulus::column-wall
                    // is zero, as assumed by default for every node
                }
                { // set the verticle flow in cement
                    axes1_as_face_normal.col(2ll) = cement_verticle_flow_vals;
                }
#pragma endregion
                return ReservoirFlowField{
                    axes1_as_face_normal,
                    axes2_as_face_normal};
            }

            template <typename Record_t, typename Grid2D_t>
            static auto create(
                const Logs::StepPropertyGrid &axes1_value,
                const Record_t &history_record,
                const Grid2D_t &grid2D)
            {
                assert(std::isnormal(history_record.rate));
                return ReservoirFlowField{
                    Logs::ZFlowRateLogFactory::create(history_record.rate, grid2D.second_coord()),
                    axes1_value,
                    grid2D};
            }
            template <typename IsPermeable_t, typename Grid2D_t>
            static auto zero_flow(
                const IsPermeable_t &is_permeable,
                const Grid2D_t &grid2D)
            {
                return horizontal_flow(0.0, is_permeable, grid2D);
            }

            template <typename IsPermeable_t, typename Grid2D_t>
            static auto horizontal_flow(
                RealType rate,
                const IsPermeable_t &is_permeable,
                const Grid2D_t &grid2D)
            {
                const Logs::StepProperty rfp{(is_permeable.log_vals * StepPropertyContainer::Constant(grid2D.first_coord().mesh_size(), rate)).eval()};
                return ReservoirFlowField{
                    Logs::ZFlowRateLogFactory::create(0.0, grid2D.second_coord()),
                    Logs::RFP_weights{Logs::StepPropertyGrid{rfp, is_permeable.grid},
                              is_permeable},
                    grid2D};
            }
        };

        struct HeatFlowField
            : public ReservoirFlowField
        {
            HeatFlowField(
                ReservoirFlowField &&flow_field,
                const RealType volumetric_heat_capacity)
                : ReservoirFlowField{std::move(flow_field)}
            {
                multiply(*this, volumetric_heat_capacity);

                axes2_as_face_normal_pos = (axes2_as_face_normal + axes2_as_face_normal.abs()) / 2.0;
                axes2_as_face_normal_neg = (axes2_as_face_normal - axes2_as_face_normal.abs()) / 2.0;

                axes1_as_face_normal_pos = (axes1_as_face_normal + axes1_as_face_normal.abs()) / 2.0;
                axes1_as_face_normal_neg = (axes1_as_face_normal - axes1_as_face_normal.abs()) / 2.0;
                {
                    assert(axes2_as_face_normal_pos.cols() == axes2_as_face_normal.cols());
                    // assert(std::all_of(axes2_as_face_normal_pos.cbegin(), axes2_as_face_normal_pos.cend(), [](const RealType v)
                    //                    { return v >= 0.0; }));
                    assert(axes2_as_face_normal_neg.cols() == axes2_as_face_normal.cols());
                    // assert(std::all_of(axes2_as_face_normal_neg.cbegin(), axes2_as_face_normal_neg.cend(), [](const RealType v)
                    //                    { return v <= 0.0; }));

                    assert(axes1_as_face_normal_pos.cols() == axes1_as_face_normal.cols());
                    // assert(std::all_of(axes1_as_face_normal_pos.cbegin(), axes1_as_face_normal_pos.cend(), [](const RealType v)
                    //                    { return v >= 0.0; }));
                    assert(axes1_as_face_normal_neg.cols() == axes1_as_face_normal.cols());
                    // assert(std::all_of(axes1_as_face_normal_neg.cbegin(), axes1_as_face_normal_neg.cend(), [](const RealType v)
                    //                    { return v <= 0.0; }));
                }
            }

            // pos and neg components of flux
            FaceValuesContainer axes1_as_face_normal_pos, axes1_as_face_normal_neg;
            FaceValuesContainer axes2_as_face_normal_pos, axes2_as_face_normal_neg;
        };

    } // Properties
} // GPN
