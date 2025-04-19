#pragma once

#include <vector>
#include <memory>

#include <Eigen/SparseCore>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/Coefficients.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace SplittingMethod
        {
            struct TemporalTerm
            {
                using Factor = Properties::Fields::Conductivity_f;

                template <typename Grid_t>
                TemporalTerm(
                    std::shared_ptr<Properties::Fields> properties,
                    const Grid_t &grid)
                    : factor{properties->capacity_vol_f} // volumes are taken into account
                {
                }

                auto DivideByTemporalStep(RealType tau) const
                {
                    return (1.0 / tau) * factor;
                }

            protected:
                Factor factor;
            };

            // single sparse matrix
            using SpMatrix = Eigen::SparseMatrix<RealType>;
            // container of such sparse matricies
            using VectSpMatrix = std::vector<SpMatrix>;

            template <typename LaplaceFactor_t>
            struct BaseSplit
            {
                BaseSplit(
                    const LaplaceFactor_t &laplace_factor,
                    size_t serial_nodes, // number of matricies
                    size_t matrix_size)  // nmbr of unknowns
                    : laplace_factor{laplace_factor},
                      its_LaplaceTerm(
                          serial_nodes,                 // number of matricies
                          Eigen::SparseMatrix<RealType>{// ctor per matrix
                                                        ptrdiff_t(matrix_size),
                                                        ptrdiff_t(matrix_size)})
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