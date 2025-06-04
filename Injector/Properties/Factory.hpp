#pragma once

#include <vector>

namespace GPN
{
    namespace Logs
    {
        struct RawDataFactory
        {
    //        using Grid_t = Grids::AxesGrid<CoordinateTypes::Z>;
#pragma region HYDRODYNAMIC-LOGS

            static auto generate_is_permeable_const(
                RealType val,
                const auto &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size, val);
                return vals;
            }

            static auto generate_is_permeable(
                const auto &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 0.0 : 1.0;
                return vals;
            }

            static auto generate_porosity(
                const auto &dual_stencils,
                const auto &is_permeable)
            {
                auto size{dual_stencils.size() - 1};
                LogValuesContainer vals(size);

                for (auto id{size - size}; id < size; ++id)
                {
                    vals[id] = (id % 2 == 1) ? 0.2 : 0.5;
                    vals[id] *= is_permeable[id];
                }
                return vals;
            }

            static auto generate_permeability(
                const auto &dual_stencils,
                const auto &is_permeable)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                {
                    vals[id] = (id % 2 == 1) ? 500 : 300;
                    vals[id] *= is_permeable[id];
                }
                return vals;
            }

            static auto generate_ext_pressure(
                const auto &dual_stencils,
                const auto &is_permeable)
            {
                constexpr RealType BarToPa = 1e5;
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                {
                    vals[id] = (id % 2 == 1) ? 130 : 150;
                    vals[id] *= is_permeable[id] * BarToPa;
                }
                return vals;
            }

            static auto generate_skin(
                const auto &dual_stencils,
                const auto &is_permeable)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                {
                    vals[id] = 0.0;
                    vals[id] *= is_permeable[id];
                }
                return vals;
            }
#pragma endregion
#pragma region HEAT-LOGS
            template <typename Container_t>
            static auto generate_solid_density(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                LogValuesContainer vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 800.0 : 1000.0;
                return vals;
            }

            template <typename Container_t>
            static auto generate_solid_specific_heatcapacity(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                LogValuesContainer vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 0.7 : 0.9;
                return vals;
            }

            template <typename Container_t>
            static auto generate_conductivity(
                const Container_t &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 200 : 1000;
                return vals;
            }
#pragma endregion

            static auto generate_rates(
                const auto &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 50.0 : 0.0;
                return vals;
            }
            
            static auto generate_temperatures(
                const auto &dual_stencils)
            {
                auto size{dual_stencils.size() - 1};
                std::vector<RealType> vals(size);

                for (auto id{size - size}; id < size; ++id)
                    vals[id] = (id % 2 == 1) ? 293.0 : 273.0;
                return vals;
            }
        };
    }
}