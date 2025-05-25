#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Logs.hpp>

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
                const Grid_t &grid)
            {
                const auto temp{make_rates(well_rate, grid)};
                return ZFlowRateLog{temp};
            }

        private:
            template <typename Grid_t>
            static StepPropertyGrid make_rates(
                RealType well_rate,
                const Grid_t &grid)
            {
                StepPropertyContainer verticle_rates_vals{
                    StepPropertyContainer::Zero(grid.mesh_size())};
                verticle_rates_vals(0ll) = well_rate;

                return {verticle_rates_vals,
                        grid};
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
                const Grid2D_t &grid)
            {
                assert(grid.second_coord.mesh_size() == log.size());
                FaceValuesContainer face_vals(grid.first_coord.dual_size(), log.size());
                face_vals.rowwise() = log.log_vals.transpose();
                return face_vals;
            }

            template <typename Grid2D_t>
            static auto flow_in_dir2(
                const Logs::StepPropertyGrid &log,
                const Grid2D_t &grid)
            {
                assert(grid.first_coord.mesh_size() == log.size());
                FaceValuesContainer face_vals(log.size(), grid.second_coord.dual_size());
                face_vals.colwise() = log.log_vals;
                return face_vals;
            }
        };

        struct ReservoirFlowField
        {
            template <typename Grid2D_t>
            ReservoirFlowField(
                const FaceValuesContainer &axes1_value,
                const FaceValuesContainer &axes2_value,
                const Grid2D_t &grid)
                : axes1_as_face_normal{axes1_value},
                  axes2_as_face_normal{axes2_value}
            {
            }

            template <typename Well_t, typename Grid2D_t>
            ReservoirFlowField(
                RealType well_rate,
                const Well_t &well,
                const Grid2D_t &grid)
                : ReservoirFlowField{
                      Logs::ZFlowRateLog{well_rate, grid.second_coord},
                      Logs::RFP{well_rate, well},
                      grid}
            {
            }

            template <typename Grid2D_t>
            ReservoirFlowField(
                const Logs::StepPropertyGrid &axes1_value,
                const Logs::StepPropertyGrid &axes2_value,
                const Grid2D_t &grid)
                : axes1_as_face_normal{FlowFieldFactory::flow_in_dir1(axes1_value, grid)},
                  axes2_as_face_normal{FlowFieldFactory::flow_in_dir2(axes2_value, grid)}
            {
            }

            template <typename Grid2D_t>
            ReservoirFlowField(
                const Logs::StepPropertyGrid &axes1_value,
                RealType well_rate,
                const Grid2D_t &grid)
                : ReservoirFlowField{
                      Logs::ZFlowRateLog{well_rate, grid.first_coord},
                      axes1_value,
                      grid}
            {
            }

            FaceValuesContainer axes1_as_face_normal;
            FaceValuesContainer axes2_as_face_normal;
        };

        struct FlowFactory
        {
            template <typename IsPermeable_t, typename Grid2D_t>
            static auto zero_flow(
                const IsPermeable_t &is_permeable,
                const Grid2D_t &grid)
            {
                return horizontal_flow(0.0, is_permeable, grid);
            }

            template <typename IsPermeable_t, typename Grid2D_t>
            static auto horizontal_flow(
                RealType rate,
                const IsPermeable_t &is_permeable,
                const Grid2D_t &grid)
            {
                const auto rfp{StepPropertyContainer::Constant(grid.first_coord.mesh_size(), rate)};
                return ReservoirFlowField{
                    FlowFieldFactory::flow_in_dir1(
                        Logs::ZFlowRateLogFactory::create(0.0, grid.second_coord),
                        grid),
                    FlowFieldFactory::flow_in_dir2(
                        Logs::RFP{Logs::StepPropertyGrid{rfp, is_permeable.grid},
                                  is_permeable},
                        grid),
                    grid};
            }
        };
    } // Properties
} // GPN
