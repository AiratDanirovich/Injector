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
            /// @brief The split is along second direction (y),
            /// the 1D problems are solved along the first direction (x)
            template <typename Grid_t, typename LaplaceFactor_t>
            struct SplitY
                : public BaseSplit<LaplaceFactor_t>
            {
                SplitY(
                    const LaplaceFactor_t &laplace_factor,
                    const Grid_t &grid)
                    : BaseSplit<LaplaceFactor_t>{
                        laplace_factor,
                          grid.second_coord.size(), // nmbr of matricies
                          grid.first_coord.size()}  // nmbr of unknowns
                {
                    assert(laplace_factor.cols() == grid.second_coord.mesh_size());
                    assert(laplace_factor.rows() == grid.first_coord.dual_size() - 2ll);
                    FillLaplaceTerm(grid);
                }

            protected:
                void FillLaplaceTerm(const Grid_t &grid)
                {
                    const auto &x_face_factor = laplace_factor.face_vals_axes1;
                    assert(x_face_factor.cols() == static_cast<std::ptrdiff_t>(its_LaplaceTerm.size()));
                    const auto &x_face_area = grid.face_area_axes1;
                    assert(x_face_factor.cols() == x_face_area.size());

#pragma omp parallel for
                    for (ptrdiff_t m_id = 0; m_id < ptrdiff_t(its_LaplaceTerm.size()); ++m_id)
                    {
                        auto &matrix = its_LaplaceTerm[m_id];
                        std::vector<Eigen::Triplet<RealType, ptrdiff_t>> tripletList;
                        tripletList.reserve(matrix.rows() * 3ull - 2ull);

                        tripletList.emplace_back(
                            0, 0,
                            // conductivity_y_bounds only contains internal boundaries of control volumes
                            x_face_factor(0, m_id) * x_face_area(m_id));
                        tripletList.emplace_back(
                            0, 1,
                            -x_face_factor(1, m_id) * x_face_area(m_id));
                        for (std::ptrdiff_t row{1ll}; row < std::ptrdiff_t(matrix.rows()) - 1ll; ++row)
                        {
                            tripletList.emplace_back(
                                row, row - 1,
                                -x_face_factor(row - 1ll, m_id) * x_face_area(m_id));
                            tripletList.emplace_back(
                                row, row,
                                (x_face_factor(row - 1ll, m_id) +
                                 x_face_factor(row, m_id)) *
                                    x_face_area(m_id));
                            tripletList.emplace_back(
                                row, row + 1,
                                -x_face_factor(row, m_id) * x_face_area(m_id));
                        }
                        auto end{std::ptrdiff_t(matrix.rows()) - 1ll};
                        tripletList.emplace_back(
                            end, end - 1ll,
                            -x_face_factor(end - 1ll, m_id) * x_face_area(m_id));
                        tripletList.emplace_back(
                            end, end,
                            x_face_factor(end - 1ll, m_id) * x_face_area(m_id));

                        matrix.setFromTriplets(tripletList.begin(), tripletList.end());
                    }
                }
            };
        } // Splitting method
    } // EqSolver
} // GPN