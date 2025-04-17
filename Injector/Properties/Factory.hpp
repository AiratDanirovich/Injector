#pragma once

#include <vector>

#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    namespace Logs
    {
        struct Factory
        {
            using Grid_t = Grids::AxesGrid<CoordinateTypes::Z>;
            
            static auto generate_permeability_StepProperty(const Grids::GridDualStencils& dual_stencils)
            {
                std::vector<RealType> vals(dual_stencils.size());

                for(ptrdiff_t id{0}; id < dual_stencils.size(); ++id)
                    vals[id] = (id % 2 == 1) ? 500 : 300;
                return vals;
            }
            
            static auto generate_porosity_StepProperty(const Grids::GridDualStencils& dual_stencils)
            {
                std::vector<RealType> vals(dual_stencils.size());

                for(ptrdiff_t id{0}; id < dual_stencils.size(); ++id)
                    vals[id] = (id % 2 == 1) ? 0.2 : 0.5;
                return vals;
            }

            static auto generate_rates_StepProperty(const Grids::GridDualStencils& dual_stencils)
            {
                std::vector<RealType> vals(dual_stencils.size());

                for(ptrdiff_t id{0}; id < dual_stencils.size(); ++id)
                    vals[id] = (id % 2 == 1) ? 50.0 : 0.0;
                return vals;
            }
        };
    }
}