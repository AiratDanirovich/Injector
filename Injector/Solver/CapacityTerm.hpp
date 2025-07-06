#pragma once

#include <Injector/Grids/Defines.h>

#include <Eigen/Core>

namespace GPN
{    
    namespace EqSolver
    {
        struct TemporalTerm
            {
                template <typename Capacity_t, typename Grid_t>
                TemporalTerm(
                    const Capacity_t &factor,
                    const Grid_t &grid)
                    : capacity{grid.volumes() * factor.values()} // volumes are taken into account
                {
                    for (auto j{0ll}; j < this->capacity.cols(); ++j)
                        for (auto i{0ll}; i < this->capacity.rows(); ++i)
                        {
                            assert(!std::isinf(this->capacity(i, j)));
                            assert(!std::isnan(this->capacity(i, j)));
                        }
                }

                auto Divide(RealType tau) const
                {
                    assert(tau > 0.0);
                    return (static_cast<RealType>(1.0) / tau) * capacity;
                }

            protected:
                // multiplied by cell volume
                const Eigen::ArrayXX<RealType> capacity;
            };

    } // EqSolver
} // GPN