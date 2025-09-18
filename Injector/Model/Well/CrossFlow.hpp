#pragma once

#include <vector>
#include <algorithm>
#include <iterator>

#include <Injector/Grids/Defines.h>

#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    /// @brief Descriptor of a single cross flow flux,
    /// from a single hole in the column hrough the cement to a permeable layer
    struct SingleCrossFlow
    {
        SingleCrossFlow(
            const Logs::StepPropertyGrid &some_log,
            const std::ptrdiff_t from_id,
            const std::ptrdiff_t to_id,
            const RealType flux)
            : // dir = 1.0 if the cross flow is directed downwards,
              // dir = -1.0 if the crossflow is directed upwards
              verticle_flux{
                  set_verticle_flux(
                      some_log, from_id, to_id,
                      flux * (from_id < to_id ? 1.0 : -1.0))},
              from_id{from_id},
              to_id{to_id}
        {
        }

        // the flux values take into account the direction of flow
        // const RealType flux;
        const std::ptrdiff_t from_id, to_id;
        const StepPropertyContainer verticle_flux;

    private:
        static StepPropertyContainer set_verticle_flux(
            const Logs::StepPropertyGrid &some_log,
            const std::ptrdiff_t from_id,
            const std::ptrdiff_t to_id,
            const RealType directed_flux)
        {
            assert((from_id >= 0ll) && (from_id < some_log.log_vals.rows()));
            assert((to_id >= 0ll) && (to_id < some_log.log_vals.rows()));
            StepPropertyContainer out{StepPropertyContainer::Zero(some_log.log_vals.rows())};
            out.middleRows(from_id, std::abs(from_id - to_id)) = directed_flux;
            return out;
        }
    };

    struct CrossFlows
    {
        CrossFlows(
            const RealType total_rate,
            const Logs::RateWeights &rate_weights,
            const std::vector<RealType> &from_coords,
            const std::vector<ptrdiff_t> &to_layers)
        //    : cross_flow_data{}
        {
        }

    protected:
        std::vector<SingleCrossFlow> cross_flow_data;

    private:
        static std::vector<SingleCrossFlow> set_cross_flow_data(
            const RealType total_rate,
            const Logs::RateWeights &rate_weights,
            const std::vector<RealType> &from_coords,
            const std::vector<ptrdiff_t> &to_layers)
        {
            assert(from_coords.size() == to_layers.size());

            const auto& grid{rate_weights.grid};

            std::vector<SingleCrossFlow> out;
            out.reserve(from_coords.size());

            for(auto i{0ull}; i < from_coords.size(); ++i)
            {

            }

        }

        static std::ptrdiff_t get_cell_id(
            const Logs::StepPropertyGrid::Grid_t &grid,
            const RealType from_coord)
        {
            const auto& dual_nodes{grid.dual_nodes};

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
            assert((dual_nodes(idx) < from_coord) && (dual_nodes(idx+1ll) > from_coord));
            
            const auto& mesh_nodes{grid.mesh_nodes};
            assert(std::abs(mesh_nodes(idx) - from_coord) < (dual_nodes(idx+1ll) - dual_nodes(idx))/2.0);

            return idx;
        }
    };

} // GPN