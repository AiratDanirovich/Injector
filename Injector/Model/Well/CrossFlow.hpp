#pragma once

#include <vector>
#include <algorithm>
#include <iterator>

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/LogsFactory.hpp>

namespace GPN
{
    namespace CrossFlow
    {
        /// @brief Descriptor of a single cross flow flux,
        /// from a single hole in the column vertically along the cement towards a permeable layer
        struct SingleCrossFlow
        {
            /// @brief Container for verticle flux in the cement due to a singel hole in the column
            /// @param some_log just some log to get the grid pointer and log_vals.size
            /// @param from_id id of cell with a hole in the column
            /// @param to_id id of cell with a hole in the sandface
            /// @param flux_amount An element of RFP profile
            SingleCrossFlow(
                const Logs::StepPropertyGrid &some_log,
                const std::ptrdiff_t from_id,
                const std::ptrdiff_t to_id,
                const RealType flux_amount = 0.0) // RFP value
                : verticle_flux{
                      StepPropertyContainer::Zero(
                          some_log.grid.dual_nodes.rows())},
                  from_id{from_id}, to_id{to_id}, dir{get_dir(from_id, to_id)}
            {
                assert((from_id >= 0ll) && (from_id < some_log.log_vals.rows()));
                assert((to_id >= 0ll) && (to_id < some_log.log_vals.rows()));

                set_flux(flux_amount);
            }

            // the flux values take into account the direction of flow
            const std::ptrdiff_t from_id, to_id;
            const RealType dir;
            const StepPropertyContainer verticle_flux;

            void set_flux(const Logs::RFP &flux)
            {
                set_flux(flux(to_id));
            }

            void set_flux(const RealType flux_amount)
            {
                const RealType directed_flux{flux_amount * dir};
                const_cast<StepPropertyContainer &>(verticle_flux).middleRows(std::min(from_id, to_id) + 1ll, std::abs(from_id - to_id)) =
                    directed_flux;
            }

#pragma region PRIVATE-METHODS
        private:
            static RealType get_dir(
                const std::ptrdiff_t from_id,
                const std::ptrdiff_t to_id)
            {
                // dir = 1.0 if the cross flow is directed downwards,
                // dir = -1.0 if the crossflow is directed upwards
                return (from_id < to_id ? 1.0 : -1.0);
            }
#pragma endregion
        };

        struct CrossFlowsHandler
        {
            CrossFlowsHandler(
                const std::vector<RealType> &from_coords,
                const std::vector<ptrdiff_t> &to_layers,
                const Logs::StepPropertyGrid &some_log)
                : cross_flow_data{
                      set_cross_flow_data(
                          some_log, from_coords, to_layers)}
            {
                // verify that the holes in the column
                // are separated by a sufficient space, e.g. 1 m
                auto this_coord{from_coords};
                std::sort(this_coord.begin(), this_coord.end());

                constexpr RealType min_spacing{1.0};
                for (auto i{1ull}; i < this_coord.size(); ++i)
                    assert(this_coord[i] >= this_coord[i - 1ull] + min_spacing);
            }

            SingleCrossFlow &operator[](const size_t i)
            {
                return const_cast<std::vector<SingleCrossFlow> &>(cross_flow_data)[i];
            }
            const SingleCrossFlow &operator[](const size_t i) const
            {
                return cross_flow_data[i];
            }

            const size_t size() const
            {
                return cross_flow_data.size();
            }

            const std::vector<SingleCrossFlow> cross_flow_data;
#pragma region PRIVATE-METHODS
        private:
            /// @brief Initialize every cross-flow with a zero-value flow
            /// @param some_log just some log to get the grid
            /// @param from_coords
            /// @param to_layers
            /// @return
            static std::vector<SingleCrossFlow> set_cross_flow_data(
                const Logs::StepPropertyGrid &some_log,
                const std::vector<RealType> &from_coords,
                const std::vector<ptrdiff_t> &to_layers)
            {
                assert(from_coords.size() == to_layers.size());

                const auto &grid{some_log.grid};

                std::vector<SingleCrossFlow> out;
                out.reserve(from_coords.size());

                for (auto i{0ull}; i < from_coords.size(); ++i)
                {
                    const auto from_cell{get_mesh_cell_id(grid, from_coords[i])};
                    const auto to_cell{get_to_cell_id(grid, to_layers[i])};
                    out.emplace_back(some_log, from_cell, to_cell, 0.0);
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

        /// @brief Assumption: multiple flows can start from the same well node,
        /// but they have to come to different layers of the reservoir.
        /// This allows to assemble the WFP from the RFP
        /// @param is_permeable
        /// @param is_perforated
        /// @param RFP_weights
        /// @param cross_flows
        /// @return
        Logs::WFP create_WFP_weights(
            const Logs::IsPerforated &is_perforated,
            const Logs::RFP &RFP_w,
            const CrossFlowsHandler &cross_flow_handler)
        {
            // first, account for true perforations
            StepPropertyContainer wfp_step_prop_grid{(is_perforated * RFP_w).log_vals};

            for (const auto &cf : cross_flow_handler.cross_flow_data)
                wfp_step_prop_grid(cf.from_id) += RFP_w(cf.to_id);

            assert(wfp_step_prop_grid.sum() == RFP_w.log_vals.sum());

            assert(is_perforated.size() == wfp_step_prop_grid.size());
            
            const auto is_damaged{Logs::IsDamagedFactory::create(
                cross_flow_handler,
                is_perforated.log_vals,
                is_perforated.grid)};

            // assert flow against perforated and damaged cells
            for (auto i{0ll}; i < wfp_step_prop_grid.size(); ++i)
            {
                assert(
                    ((is_perforated(i) != is_damaged(i)) //&&
                    // (wfp_step_prop_grid(i) != 0.0)
                    ) ||
                    ((is_perforated(i) == 0.0) && (is_damaged(i) == 0.0) &&
                     (wfp_step_prop_grid(i) == 0.0)));
            }

            return Logs::WFPFactory::create_from_container<Logs::WFP>(
                wfp_step_prop_grid, is_perforated + is_damaged);
        }

        struct CrossFlows
        {
            CrossFlows(
                const std::vector<RealType> &from_coords,
                const std::vector<ptrdiff_t> &to_layers,
                const Logs::RFP& rfp,
                const Logs::IsPerforated &is_perforated)
                : cross_flow_handler{
                      from_coords, to_layers, is_perforated},
                  is_perforated{is_perforated}
            {
                set_flux(rfp);
            }
            CrossFlows(
                const std::vector<RealType> &from_coords,
                const std::vector<ptrdiff_t> &to_layers,
                const Logs::IsPerforated &is_perforated)
                : cross_flow_handler{
                      from_coords, to_layers, is_perforated},
                  is_perforated{is_perforated}
            {
            }

            void set_flux(const Logs::RFP &RFP_w)
            {
                const_cast<StepPropertyContainer&>(rfp) = RFP_w.log_vals;
                for (auto id{0ull}; id < cross_flow_handler.size(); ++id)
                    const_cast<CrossFlowsHandler &>(cross_flow_handler)[id].set_flux(RFP_w);

                const_cast<StepPropertyContainer &>(verticle_flux_in_cement) = get_cement_verticle_flux(RFP_w);

                const_cast<StepPropertyContainer&>(wfp) = 
                    create_WFP_weights(
                        is_perforated,
                        RFP_w,
                        cross_flow_handler).log_vals;
                
                const_cast<StepPropertyContainer&>(verticle_flux_in_well) = get_well_verticle_flux(wfp);
            }

            const StepPropertyContainer get_verticle_flux() const
            {
                return verticle_flux_in_cement;
            }

            const CrossFlowsHandler cross_flow_handler;
            const StepPropertyContainer rfp, wfp;
            const StepPropertyContainer verticle_flux_in_cement, verticle_flux_in_well;

#pragma region PRIVET-METHODS
        private:
            const StepPropertyContainer get_cement_verticle_flux(
                const Logs::RFP &RFP_w) const
            {
                StepPropertyContainer out{StepPropertyContainer::Zero(RFP_w.grid.dual_nodes.rows())};
                for (auto i{0ull}; i < cross_flow_handler.size(); ++i)
                    out += cross_flow_handler[i].verticle_flux;
                return out;
            }
            
            const StepPropertyContainer get_well_verticle_flux(
                StepPropertyContainer wfp) const
            {
                const RealType rate{wfp.sum()};
                StepPropertyContainer out(StepPropertyContainer::Zero(wfp.rows() + 1ll));
                std::partial_sum(wfp.cbegin(), wfp.cend(), out.begin() + 1ll, std::plus<RealType>{});
                return rate - out;
            }

            const Logs::IsPerforated &is_perforated;
#pragma endregion
        };
    } // CrossFlow
} // GPN