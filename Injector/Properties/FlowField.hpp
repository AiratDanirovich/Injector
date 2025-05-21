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
            template <typename Grid_t>
            ZFlowRateLog(
                RealType well_rate,
                const Grid_t &grid)
                : StepPropertyGrid{make_rates(well_rate, grid)},
                  AssertNonNegative{make_rates(well_rate, grid)}
            {
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

    namespace Properties
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
            template <typename Well_t, typename Grid2D_t>
            ReservoirFlowField(
                RealType well_rate,
                const Well_t &well,
                const Grid2D_t &grid)
                : ReservoirFlowField{
                      Logs::ZFlowRateLog{well_rate, grid.second_coord},
                      Logs::RFP{well_rate,
                                well,
                                grid.first_coord},
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
    } // Properties
} // GPN
