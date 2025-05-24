#pragma once

#include <vector>
#include <algorithm>
#include <numeric>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids2D.hpp>

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

            static auto create_cartesian_grid_2D_ptr(ptrdiff_t n)
            {
                auto stencils{Factory::generate_dual_grid_stencils_uniform(0, 1, n)};

                auto nodes{GridDual{stencils}};

                auto x_grid{
                    AxesGrid<CoordinateTypes::X>{nodes}};

                auto y_grid{
                    AxesGrid<CoordinateTypes::Y>{nodes}};

                return std::make_shared<StructuredXYGrid2D>(x_grid, y_grid);
            }

            static auto create_cylinder_grid_2D_ptr(const Box &box, ptrdiff_t n1, ptrdiff_t n2)
            {
                auto z_stencils{Factory::generate_dual_grid_stencils_uniform(box.axes1, n1)};
                auto r_stencils{Factory::generate_dual_grid_stencils_uniform(box.axes2, n2)};
                return create_cylinder_grid_2D_ptr(z_stencils, r_stencils);
            }

            template <typename CoordinateType_t>
            static auto create_axes(const auto &stencils)
            {
                auto nodes{GridDual{stencils}};
                return 
                    AxesGrid<CoordinateType_t>{nodes};
            }

            static auto create_cylinder_grid_2D_ptr(const auto &z_stencils, const auto &r_stencils)
            {
                //    auto z_nodes{GridDual{z_stencils}};
                auto z_grid{create_axes<CoordinateTypes::Z>(z_stencils)};
                //        AxesGrid<CoordinateTypes::Z>{z_nodes}};

                auto r_grid{create_axes<CoordinateTypes::R_CylCoord>(r_stencils)};
                // auto r_nodes{GridDual{r_stencils}};
                // auto r_grid{
                //     AxesGrid<CoordinateTypes::R_CylCoord>{r_nodes}};

                return std::make_shared<StructuredCylinderGrid2DAxisymmetric>(z_grid, r_grid);
            }
        };

        struct CylinderGridFactory
        {
            CylinderGridFactory(const Box &box, ptrdiff_t n1, ptrdiff_t n2)
                : CylinderGridFactory(
                      Factory::generate_dual_grid_stencils_uniform(box.axes1, n1),
                      Factory::generate_dual_grid_stencils_uniform(box.axes2, n2))
            {
            }

            CylinderGridFactory(const auto &z_stencils, const auto &r_stencils)
                : grid2D{
                      Grids::Factory::create_cylinder_grid_2D_ptr(
                          z_stencils, r_stencils)}
            {
            }

            const auto grid() const
            {
                return grid2D;
            }

        protected:
            cptr<Grids::StructuredCylinderGrid2DAxisymmetric> grid2D;
        };
    }
}