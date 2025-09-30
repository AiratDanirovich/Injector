#pragma once

#include <vector>
#include <algorithm>
#include <iterator>

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    namespace CrossFlow
    {
        /// @brief Descriptor of a single cross flow flux,
        /// from a single hole in the column hrough the cement to a permeable layer
        struct SingleCrossFlow
        {
            /// @brief
            /// @param some_log
            /// @param from_id
            /// @param to_id
            /// @param flux_fraction An element of normalized RFP profile
            SingleCrossFlow(
                const Logs::StepPropertyGrid &some_log,
                const std::ptrdiff_t from_id,
                const std::ptrdiff_t to_id,
                const RealType flux_fraction)
                : verticle_flux{
                      set_verticle_flux(
                          some_log, from_id, to_id,
                          flux_fraction)},
                  from_id{from_id}, to_id{to_id}
            {
            }

            // the flux values take into account the direction of flow
            // const RealType flux;
            const std::ptrdiff_t from_id, to_id;
            const StepPropertyContainer verticle_flux;
#pragma region PRIVATE-METHODS
        private:
            static StepPropertyContainer set_verticle_flux(
                const Logs::StepPropertyGrid &some_log,
                const std::ptrdiff_t from_id,
                const std::ptrdiff_t to_id,
                const RealType flux_fraction)
            {
                // dir = 1.0 if the cross flow is directed downwards,
                // dir = -1.0 if the crossflow is directed upwards
                const RealType dir{(from_id < to_id ? 1.0 : -1.0)};
                const RealType directed_normalized_flux{flux_fraction * dir};

                assert((from_id >= 0ll) && (from_id < some_log.log_vals.rows()));
                assert((to_id >= 0ll) && (to_id < some_log.log_vals.rows()));
                StepPropertyContainer out{StepPropertyContainer::Zero(some_log.grid.dual_nodes.rows())};
                out.middleRows(std::min(from_id, to_id) + 1ll, std::abs(from_id - to_id)) = directed_normalized_flux;
                return out;
            }
#pragma endregion
        };

        struct CrossFlows
        {
            CrossFlows(
                const Logs::RFP &RFP_weights,
                const std::vector<RealType> &from_coords,
                const std::vector<ptrdiff_t> &to_layers)
                : cross_flow_data{
                      set_cross_flow_data(
                          RFP_weights,
                          from_coords, to_layers)}
            {
                normalized_verticle_flux = set_verticle_flux(RFP_weights);

                auto this_coord{from_coords};
                std::sort(this_coord.begin(), this_coord.end());

                constexpr RealType min_spacing{1.0};
                for (auto i{1ull}; i < this_coord.size(); ++i)
                    assert(this_coord[i] >= this_coord[i - 1ull] + min_spacing);
            }

            const StepPropertyContainer get_normalized_verticle_flux() const
            {
                return normalized_verticle_flux;
            }

            const StepPropertyContainer get_verticle_flux(
                const RealType total_flux) const
            {
                return total_flux * get_normalized_verticle_flux();
            }

            const std::vector<SingleCrossFlow> cross_flow_data;
#pragma region PRIVETA-METHODS
        private:
            StepPropertyContainer normalized_verticle_flux;
            const StepPropertyContainer set_verticle_flux(const Logs::RFP &RFP_weights) const
            {
                StepPropertyContainer out{StepPropertyContainer::Zero(RFP_weights.grid.dual_nodes.rows())};
                for (auto i{0ull}; i < cross_flow_data.size(); ++i)
                    out += cross_flow_data[i].verticle_flux;
                return out;
            }

            static std::vector<SingleCrossFlow> set_cross_flow_data(
                const Logs::RFP &RFP_weights,
                const std::vector<RealType> &from_coords,
                const std::vector<ptrdiff_t> &to_layers)
            {
                assert(from_coords.size() == to_layers.size());

                const auto &grid{RFP_weights.grid};

                std::vector<SingleCrossFlow> out;
                out.reserve(from_coords.size());

                for (auto i{0ull}; i < from_coords.size(); ++i)
                {
                    const auto from_cell{get_mesh_cell_id(grid, from_coords[i])};
                    const auto to_cell{get_to_cell_id(grid, to_layers[i])};
                    out.emplace_back(RFP_weights, from_cell, to_cell, RFP_weights(to_cell));
                }

                return out;
            }

            static std::ptrdiff_t get_mesh_cell_id(
                const Logs::StepPropertyGrid::Grid_t &grid,
                const RealType from_coord)
            {
                const auto &dual_nodes{grid.dual_nodes};

                assert(from_coord > dual_nodes(0ll));
                assert(from_coord < dual_nodes.tail<1ll>()(0ll));

                const auto it{
                    std::upper_bound(
                        dual_nodes.begin(),
                        dual_nodes.end(),
                        from_coord)};
                // idx in mesh_nodes container of the cell including the from_coord
                const auto idx{std::distance(dual_nodes.begin(), it) - 1ll};

                assert(idx >= 0ll);
                assert((dual_nodes(idx) < from_coord) && (dual_nodes(idx + 1ll) > from_coord));

                const auto &mesh_nodes{grid.mesh_nodes};
                assert(std::abs(mesh_nodes(idx) - from_coord) < (dual_nodes(idx + 1ll) - dual_nodes(idx)) / 2.0);

                return idx;
            }

            static std::ptrdiff_t get_to_cell_id(
                const Logs::StepPropertyGrid::Grid_t &grid,
                const std::ptrdiff_t to_id)
            {
                const auto &dual_stencils{grid.dual_stencils};

                // middle coordinate of the permeable layer that
                // accepts the cross flow
                const RealType to_coord{(dual_stencils(to_id) + dual_stencils(to_id + 1ll)) / 2.0};

                return get_mesh_cell_id(grid, to_coord);
            }
#pragma endregion
        };
    } // CrossFlow
} // GPN