#pragma once

#include <Eigen/Core>
// #include <Eigen/Dense>

#include <Injector/Grids/Defines.h>

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

                template <typename StructuredGrid2D_t>
                static auto FillWithZeros(
                    const StructuredGrid2D_t &grid)
                {
                    return FillWithConst(grid, (RealType)0.0);
                }

                template <typename StructuredGrid2D_t>
                static auto FillWithConst(
                    const StructuredGrid2D_t &grid,
                    RealType val)
                {
                    State_Container cur_state{
                        grid.first_coord.size(),
                        grid.second_coord.size()};
                    cur_state.fill(val);

                    return State2D{cur_state};
                }

                template <
                    typename StructuredGrid2D_t,
                    typename Functor>
                static auto FillWithFunctor(
                    const StructuredGrid2D_t &grid,
                    const Functor &f, RealType initial_moment = 0.0)
                {
                    State_Container cur_state{
                        grid.first_coord.size(),
                        grid.second_coord.size()};

                    for (std::ptrdiff_t j = 0; j < cur_state.outerSize(); ++j)
                        for (std::ptrdiff_t i = 0; i < cur_state.innerSize(); ++i)
                        {
                            cur_state(i, j) = f(
                                grid.first_coord[i],
                                grid.second_coord[j],
                                initial_moment);
                        }

                    return State2D{cur_state};
                }

                State2D(const State2D &) noexcept = default;
                State2D(State2D &&) noexcept = default;

                RealType operator()(auto i, auto j) const
                {
                    return cur_state(i, j);
                }

                RealType &operator()(auto i, auto j)
                {
                    return cur_state(i, j);
                }

            public:
                State_Container cur_state;
            };
        } // Coefficients
    } // EqSolver
} // GPN
