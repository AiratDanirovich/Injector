#pragma once

#include <cmath>
#include <Eigen/Core>

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
                        grid.first_coord.mesh_size(),
                        grid.second_coord.mesh_size()};
                    cur_state.fill(val);

                    return cur_state;
                }

                template <
                    typename StructuredGrid2D_t,
                    typename Functor>
                static auto FillWithFunctor(
                    const StructuredGrid2D_t &grid,
                    const Functor &f, RealType initial_moment = 0.0)
                {
                    State_Container cur_state{
                        grid.first_coord.mesh_size(),
                        grid.second_coord.mesh_size()};

                    for (std::ptrdiff_t j = 0; j < cur_state.outerSize(); ++j)
                        for (std::ptrdiff_t i = 0; i < cur_state.innerSize(); ++i)
                        {
                            cur_state(i, j) = f(
                                grid.first_coord.coordinate(i),
                                grid.second_coord.coordinate(j),
                                initial_moment);
                        }

                    for (std::ptrdiff_t j = 0; j < cur_state.outerSize(); ++j)
                        for (std::ptrdiff_t i = 0; i < cur_state.innerSize(); ++i)
                        {
                            assert(!std::isinf(cur_state(i, j)));
                            assert(!std::isnan(cur_state(i, j)));
                        }

                    return State2D{cur_state};
                }

                State2D(const State2D &) noexcept = default;
                State2D(State2D &&) noexcept = default;
                State2D() = delete;

                RealType operator()(auto i, auto j) const
                {
                    return cur_state(i, j);
                }

                RealType &operator()(auto i, auto j)
                {
                    return cur_state(i, j);
                }

                auto rows() const
                {
                    return cur_state.rows();
                }
                auto cols() const
                {
                    return cur_state.cols();
                }

                operator const State_Container &() const { return cur_state; }

            public:
                State_Container cur_state;
            };
        } // Coefficients
    } // EqSolver
} // GPN
