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
            
            template<typename Container_t>
            static auto generate_permeability_StepProperty(
                const Container_t& dual_stencils)
            {
                auto size{dual_stencils.size()-1};
                std::vector<RealType> vals(size);

                for(auto id{size-size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 500 : 300;
                return vals;
            }
            
            template<typename Container_t>
            static auto generate_porosity_StepProperty(
                const Container_t& dual_stencils)
            {
                auto size{dual_stencils.size()-1};
                std::vector<RealType> vals(size);

                for(auto id{size-size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 0.2 : 0.5;
                return vals;
            }
            
            template<typename Container_t>
            static auto generate_is_permeable_StepProperty(
                const Container_t& dual_stencils)
            {
                auto size{dual_stencils.size()-1};
                std::vector<RealType> vals(size);

                for(auto id{size-size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 0.0 : 1.0;
                return vals;
            }

            
            template<typename Container_t>
            static auto generate_conductivity_StepProperty(
                const Container_t& dual_stencils)
            {
                auto size{dual_stencils.size()-1};
                std::vector<RealType> vals(size);

                for(auto id{size-size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 200 : 1000;
                return vals;
            }

            static auto generate_rates_StepProperty(const Grids::GridDualStencils& dual_stencils)
            {
                auto size{dual_stencils.size()-1};
                std::vector<RealType> vals(size);

                for(auto id{size-size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 50.0 : 0.0;
                return vals;
            }
        };
    }
}