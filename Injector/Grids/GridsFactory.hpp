#pragma once

#include <vector>
#include <algorithm>
#include <numeric>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace Grids
    {
        struct Factory
        {
            static auto generate_dual_grid_stencils_uniform(RealType a, RealType b, ptrdiff_t n)
            {
                assert(n > 1);
                RealType step{(b - a) / (n - 1)};

                std::vector<RealType> result(n);

                result.front() = a;
                for (ptrdiff_t idx{1}; idx < n - 1; ++idx)
                    result[idx] = result[idx - 1] + step;
                result.back() = b;

                return result;
            }

            static auto generate_dual_grid_steps(const std::vector<RealType> &data)
            {
                assert(data.size() > 1ull);
                std::vector<RealType> out(data.size() - 1ull);
                for (auto i{1ull}; i < data.size(); ++i)
                    out[i - 1] = data[i] - data[i - 1];
                return out;
            }

            static auto generate_dual_grid_stencils_from_steps(RealType zTop, std::vector<RealType> thickness)
            {
                thickness.insert(thickness.begin(), zTop);
                std::vector<RealType> z_stencils(thickness.size(), 0.0);
                std::partial_sum(thickness.begin(), thickness.end(), z_stencils.begin(), std::plus<RealType>());

                return z_stencils;
            }

            static auto generate_dual_grid_stencils_uniform(const Segment &axes, ptrdiff_t n)
            {
                return generate_dual_grid_stencils_uniform(axes.start, axes.end, n);
            }

            template <typename CoordinateType_t>
            static auto create_axes(const auto &stencils)
            {
                auto nodes{GridDual{stencils}};
                return AxesGrid<CoordinateType_t>{nodes};
            }

        };
    }
}