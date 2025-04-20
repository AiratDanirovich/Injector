#pragma once

#include <vector>
#include <memory>

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include <Injector/Grids/Defines.h>
// #include <Injector/Properties/Coefficients.hpp>
#include <Injector/Solver/State2D.hpp>

#include <Injector/Solver/SplittingMethod/BaseSplit.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace SplittingMethod
        {
            /// @brief The split is along the first direction (x),
            /// the 1D problems are solved along the second direction (y)
            template <typename Grid_t, typename LaplaceFactor_t>
            struct SplitX
                : public BaseSplit<LaplaceFactor_t>
            {
                using BaseSplit<LaplaceFactor_t>::laplace_factor;
                using BaseSplit<LaplaceFactor_t>::LaplaceTerm;
                using BaseSplit<LaplaceFactor_t>::size;

                SplitX(
                    const LaplaceFactor_t &laplace_factor,
                    const Grid_t &grid)
                    : BaseSplit<LaplaceFactor_t>{
                          laplace_factor,
                          grid.first_coord.size(),  // nmbr of matricies
                          grid.second_coord.size()} // nmbr of unknowns
                {
                    assert(laplace_factor.rows() == grid.first_coord.mesh_size());
                    assert(laplace_factor.cols() == grid.second_coord.dual_size() - 2ll);
                    FillLaplaceTerm(grid);
                }

            protected:
                void FillLaplaceTerm(const Grid_t &grid)
                {
                    const auto &y_face_factor = laplace_factor.face_vals_axes2;
                    assert(y_face_factor.rows() == static_cast<std::ptrdiff_t>(size()));
                    const auto &y_face_area = grid.face_area_axes2;
                    assert(y_face_factor.rows() == y_face_area.size());

#pragma omp parallel for
                    // m_id --- matrix id or stripe id. Problem is solved along every stripe
                    // m_id == fixed row in 2D grid ArrayXX representation
                    for (std::ptrdiff_t m_id = 0;
                         m_id < static_cast<std::ptrdiff_t>(size());
                         ++m_id)
                    {
                        // 3-diag matrix, rows() x rows() size
                        auto &matrix = LaplaceTerm(m_id);
                        std::vector<Eigen::Triplet<RealType, ptrdiff_t>> tripletList;
                        tripletList.reserve(matrix.rows() * 3ull - 2ull);

                        tripletList.emplace_back(
                            0, 0,
                            y_face_factor(m_id, 0) * y_face_area(m_id));
                        tripletList.emplace_back(
                            0, 1,
                            -y_face_factor(m_id, 1) * y_face_area(m_id));
                        for (ptrdiff_t row = 1; row < ptrdiff_t(matrix.rows()) - 1; ++row)
                        {
                            tripletList.emplace_back(
                                row, row - 1,
                                -y_face_factor(m_id, row - 1) * y_face_area(m_id));
                            tripletList.emplace_back(
                                row, row,
                                (y_face_factor(m_id, row - 1) +
                                 y_face_factor(m_id, row)) *
                                    y_face_area(m_id));
                            tripletList.emplace_back(
                                row, row + 1,
                                -y_face_factor(m_id, row) * y_face_area(m_id));
                        }
                        auto end{std::ptrdiff_t(matrix.rows()) - 1ll};

                        tripletList.emplace_back(
                            end, end - 1ll,
                            -y_face_factor(m_id, end - 1ll) * y_face_area(m_id));
                        tripletList.emplace_back(
                            end, end,
                            y_face_factor(m_id, end - 1ll) * y_face_area(m_id));

                        matrix.setFromTriplets(tripletList.begin(), tripletList.end());
                    }
                }
            };
        } // SplittingMethod
    } // EqSolver
} // GPN
