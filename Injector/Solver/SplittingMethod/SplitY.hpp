#pragma once

#include <vector>
#include <memory>

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Coefficients.hpp>
#include <Injector/Solver/State2D.hpp>

#include <Injector/Solver/SplittingMethod/BaseSplit.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace SplittingMethod
        {
            struct SplitY
                : public BaseSplit
            {
                template <typename Grid_t>
                SplitY(
                    std::shared_ptr<Properties::Fields> properties,
                    const Grid_t &grid)
                    : BaseSplit{
                        properties, 
                        grid.second_coord.size(), 
                        grid.first_coord.size()}
                {
                    FillLaplaceTerm(grid);
                }

            protected:
                template <typename Grid_t>
                void FillLaplaceTerm(const Grid_t &grid)
                {
                    Conductivity_f conductivity_x_bounds{
                        set_conductivity_x_bounds(properties->conductivity_f, grid)};

                    const auto &gr_x = grid.first_coord;
                    const auto &y_vol = grid.second_coord.volumes();

#pragma omp parallel for
                    for (ptrdiff_t m_id = 0; m_id < ptrdiff_t(its_LaplaceTerm.size()); ++m_id)
                    {
                        auto &matrix = its_LaplaceTerm[m_id];
                        std::vector<Eigen::Triplet<RealType, ptrdiff_t>> tripletList;
                        tripletList.reserve(matrix.rows() * 3ull - 2ull);

                        tripletList.emplace_back(
                            0, 0,
                            // conductivity_y_bounds only contains internal boundaries of control volumes
                            conductivity_x_bounds(0, m_id) / gr_x.step(0) * y_vol(0));
                        tripletList.emplace_back(
                            0, 1,
                            -conductivity_x_bounds(1, m_id) / gr_x.step(0) * y_vol(0));
                        for (ptrdiff_t row = 1; row < ptrdiff_t(matrix.rows()) - 1; ++row)
                        {
                            tripletList.emplace_back(
                                row, row - 1,
                                -conductivity_x_bounds(row-1, m_id) / gr_x.step(row - 1) * y_vol(row));
                            tripletList.emplace_back(
                                row, row,
                                (conductivity_x_bounds(row-1, m_id) / gr_x.step(row - 1) +
                                 conductivity_x_bounds(row, m_id) / gr_x.step(row)) *
                                    y_vol(row));
                            tripletList.emplace_back(
                                row, row + 1,
                                -conductivity_x_bounds(row, m_id) / gr_x.step(row) * y_vol(row));
                        }
                        ptrdiff_t end = ptrdiff_t(matrix.rows()) - 1;
                        tripletList.emplace_back(
                            end, end - 1,
                            -conductivity_x_bounds(end-1, m_id) / gr_x.step(end-1) * y_vol(end));
                        tripletList.emplace_back(
                            end, end,
                            conductivity_x_bounds(end - 1, m_id) / gr_x.step(end-1) * y_vol(end));

                        matrix.setFromTriplets(tripletList.begin(), tripletList.end());
                    }
                }

                /**
                 * @brief Average conductivity in x-direction, to get values on x-boundaries
                 * of control volumes. Each volume has two boundaries, 
                 * but two boundaries at the domain boundary are excluded, 
                 * and the number of boundaries in container is
                 * by one less than the number of control volumes. 
                 *
                 * @tparam Grid_t
                 * @param conductivity_nodes
                 * @param grid
                 * @return auto
                 */
                template <typename Grid_t>
                Conductivity_f set_conductivity_x_bounds(
                    const Conductivity_f &conductivity_nodes,
                    const Grid_t &grid) const
                {
                    Conductivity_f conductivity_x_bounds{
                        grid.X_nodes.size() - 1,
                        grid.Y_nodes.size()};

#pragma omp parallel for
                    for (ptrdiff_t j = 0; j < ptrdiff_t(conductivity_x_bounds.cols()); ++j)
                        for (ptrdiff_t i = 1; i < ptrdiff_t(conductivity_x_bounds.rows()) - 1; ++i)
                            conductivity_x_bounds(i, j) =
                                (conductivity_nodes(i + 1, j) + conductivity_nodes(i, j)) / 2.0;

                    return conductivity_x_bounds;
                }
            };

        } // Splitting method
    } // EqSolver
} // GPN