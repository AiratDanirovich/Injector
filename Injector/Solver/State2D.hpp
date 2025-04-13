#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
// #include <Injector/Grids/UniformGrid2D.h>

namespace GPN
{
    namespace EqSolver
    {
        namespace State
        {
            struct State2D
            {
                using State_Container =
                    Eigen::ArrayXX<RealType>;

                State2D(
                    const State_Container &cur_state)
                    : cur_state{cur_state}
                {
                }

                template<typename StructuredGrid2D_t>
                static State2D FillWithZeros(
                    const StructuredGrid2D_t &grid)
                {
                    State_Container cur_state{
                        grid.first_coord.size(),
                        grid.second_coord.size()};
                    cur_state.fill(0.0);

                    return {cur_state};
                }

                template<
                    typename StructuredGrid2D_t, 
                    typename Functor
                >
                static State2D FillWithFunctor(
                    const StructuredGrid2D_t &grid, 
                    const Functor &f, RealType initial_moment)
                {
                    State_Container cur_state{
                        grid.first_coord.size(),
                        grid.second_coord.size()};

                    for (ptrdiff_t j = 0; j < cur_state.outerSize(); ++j)
                        for (ptrdiff_t i = 0; i < cur_state.innerSize(); ++i)
                        {
                            cur_state(i, j) = f(grid.first_coord[i], grid.second_coord[j], initial_moment);
                        }

                    return {cur_state};
                }

                State2D(const State2D &) noexcept = default;
                State2D(State2D &&) noexcept = default;

                RealType operator()(ptrdiff_t i, ptrdiff_t j) const
                {
                    return cur_state(i, j);
                }

                RealType &operator()(ptrdiff_t i, ptrdiff_t j)
                {
                    return cur_state(i, j);
                }

            public:
                State_Container cur_state;
            };

        } // Coefficients
    } // EqSolver
} // GPN