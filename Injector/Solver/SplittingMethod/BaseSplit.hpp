#pragma once

#include <vector>
#include <memory>

#include <Eigen/SparseCore>

#include <Injector/Grids/Defines.h>
//#include <Injector/Properties/Coefficients.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace SplittingMethod
        {
            // single sparse matrix
            using SpMatrix = Eigen::SparseMatrix<RealType>;
            // container of such sparse matricies
            using VectSpMatrix = std::vector<SpMatrix>;

            template <typename LaplaceFactor_t>
            struct BaseSplit
            {
                BaseSplit(
                    const LaplaceFactor_t &laplace_factor,
                    std::ptrdiff_t serial_nodes, // number of matricies
                    std::ptrdiff_t matrix_size)  // nmbr of unknowns
                    : laplace_factor{laplace_factor},
                      its_LaplaceTerm(
                          serial_nodes,                 // number of matricies
                          Eigen::SparseMatrix<RealType>{// ctor per matrix
                                                        matrix_size,
                                                        matrix_size})
                {
                    for (auto &m : its_LaplaceTerm)
                        m.reserve(matrix_size * 3ull - 2ull);
                }

                const SpMatrix &LaplaceTerm(size_t i) const
                {
                    return its_LaplaceTerm[i];
                }

            protected:
                VectSpMatrix its_LaplaceTerm;
                // laplace_factor only contains internal boundaries of control volumes
                const LaplaceFactor_t &laplace_factor;
            };
        } // SplittingMethod
    } // EqSolver
} // GPN