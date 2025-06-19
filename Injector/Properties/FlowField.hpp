#pragma once
#include <numeric>
#include <algorithm>
#include <cmath>

#include <iostream>

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
        struct FlowFieldFactory // FaceInterpolatedField_1D
        {
            template <typename Grid2D_t>
            static auto flow_in_dir1(
                const Logs::StepPropertyGrid &log,
                const Grid2D_t &grid2D)
            {
                assert(grid2D.second_coord.mesh_size() == log.size());
                FaceValuesContainer face_vals(grid2D.first_coord.dual_size(), log.size());
                face_vals.rowwise() = log.log_vals.transpose();
                return face_vals;
            }

            template <typename Grid2D_t>
            static auto flow_in_dir2(
                const Logs::StepPropertyGrid &log,
                const Grid2D_t &grid2D)
            {
                assert(grid2D.first_coord.mesh_size() == log.size());
                FaceValuesContainer face_vals(log.size(), grid2D.second_coord.dual_size());
                face_vals.colwise() = log.log_vals;
                return face_vals;
            }
        };

        struct ReservoirFlowField
        {
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
#pragma region AXES2-AS-FACENORMAL
                const auto rfp{Logs::RFPFactory::create_from_well(history_record, well)};
                auto axes2_as_face_normal{FlowFieldFactory::flow_in_dir2(rfp, grid2D)};
                axes2_as_face_normal.col(0ll) = 0.0;                          // boundary condition, zero flux at the axis of symmetry
                axes2_as_face_normal.col(1ll) = well.get_WFP(history_record); // flow at the tube radius
#pragma endregion
#pragma region AXES1-AS-FACENORMAL
                FaceValuesContainer axes1_as_face_normal{
                    FaceValuesContainer::Zero(
                        grid2D.first_coord.dual_size(),
                        grid2D.second_coord.mesh_size())};

                { // set the flow in tube
                    const auto wfp{Logs::WFPFactory::create_from_well(history_record, well)};
                    // ref to log vals as Eigen::ArrayX container
                    const auto &wfp_vals{wfp.log_vals};
                    // cumsum of rfp flow rates
                    LogValuesContainer cum_sum{LogValuesContainer::Zero(grid2D.first_coord.dual_size())};
                    std::partial_sum(wfp_vals.cbegin(), wfp_vals.cend(), cum_sum.begin() + 1ll, std::plus<RealType>());
                    // leftover flowrate along the well
                    const LogValuesContainer z_flow{history_record.rate - cum_sum};

                    axes1_as_face_normal.col(0ll) = z_flow;
                }
                { // set the flow in cement
                    const auto rfp{Logs::RFPFactory::create_from_well(history_record, well)};
                    // ref to log vals as Eigen::ArrayX container
                    const auto &rfp_vals{rfp.log_vals};
                    // cumsum of rfp flow rates
                    LogValuesContainer cum_sum{LogValuesContainer::Zero(grid2D.first_coord.dual_size())};
                    std::partial_sum(rfp_vals.cbegin(), rfp_vals.cbegin() + well.top_collector_cell_id, cum_sum.begin() + 1ll, std::plus<RealType>());

                    std:: cout << "cum_sum:\n"
                    << cum_sum << std::endl;
                    // leftover flowrate along the well
                    const LogValuesContainer& z_flow{cum_sum};
                    
                    axes1_as_face_normal.col(1ll) = -z_flow;
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
                    Logs::ZFlowRateLogFactory::create(history_record.rate, grid2D.second_coord),
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
                const auto rfp{(is_permeable.log_vals * StepPropertyContainer::Constant(grid2D.first_coord.mesh_size(), rate)).eval()};
                return ReservoirFlowField{
                    Logs::ZFlowRateLogFactory::create(0.0, grid2D.second_coord),
                    Logs::RFP{Logs::StepPropertyGrid{rfp, is_permeable.grid},
                              is_permeable},
                    grid2D};
            }
        };

    } // Properties
} // GPN
