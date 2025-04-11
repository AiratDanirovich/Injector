#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids.h>

namespace GPN
{
    namespace Grids
    {
        struct Factory
        {
            static auto generate_stencils_uniform(float_t a, float_t b, ptrdiff_t n)
            {
                float_t step{(b - a) / (n - 1)};

                std::vector<float_t> result(n);

                result[0] = a;
                for (ptrdiff_t idx{1}; idx < n; ++idx)
                    result[idx] = result[idx - 1] + step;

                return result;
            }
        };
    }
}