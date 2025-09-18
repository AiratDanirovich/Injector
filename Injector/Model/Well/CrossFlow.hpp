#pragma once

#include <vector>

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
              veticle_flux{
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
        CrossFlows()
        {
        }

    protected:
        std::vector<SingleCrossFlow> cross_flow_data;
    };

} // GPN