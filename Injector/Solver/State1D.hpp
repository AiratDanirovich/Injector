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
            struct State1D
            {
                using State_Container =
                    Eigen::ArrayX<RealType>;

                State1D(
                    const State_Container &cur_state)
                    : cur_state{cur_state}
                {
                }
                // State1D(
                //     State_Container &&cur_state)
                //     : cur_state{std::move(cur_state)}
                // {
                // }
                // State1D(
                //     State_Container cur_state)
                //     : cur_state{cur_state}
                // {
                // }

                template <typename StructuredGrid1D_t>
                static auto FillWithZeros(
                    const StructuredGrid1D_t &grid_z)
                {
                    return FillWithConst(grid_z, (RealType)0.0);
                }

                template <typename StructuredGrid1D_t>
                static auto FillWithConst(
                    const StructuredGrid1D_t &grid_z,
                    RealType val)
                {
                    State_Container cur_state{
                        grid_z.mesh_size()};
                    cur_state.fill(val);

                    return cur_state;
                }

                template <
                    typename StructuredGrid1D_t>
                static auto FillWithFunctor(
                    const StructuredGrid1D_t &grid2D,
                    const InitialConditions::ICFunctorBase1D &f,
                    RealType initial_moment = 0.0)
                {
                    State_Container cur_state{
                        grid2D.first_coord().mesh_size()};

                        for (std::ptrdiff_t i{0ll}; i < cur_state.innerSize(); ++i)
                        {
                            cur_state(i) = f(
                                i, //grid2D.first_coord.coordinate(i),
                                initial_moment);
                        }

                        for (std::ptrdiff_t i = 0; i < cur_state.innerSize(); ++i)
                        {
                            assert(!std::isinf(cur_state(i)));
                            assert(!std::isnan(cur_state(i)));
                        }

                    return State1D{cur_state};
                }

                State1D(const State1D &) noexcept = default;
                State1D(State1D &&) noexcept = default;
                State1D() = delete;

                RealType operator()(auto i) const
                {
                    return cur_state(i);
                }

                RealType &operator()(auto i)
                {
                    return cur_state(i);
                }

                auto rows() const
                {
                    return cur_state.rows();
                }
                constexpr auto cols() const
                {
                    return 1ll;
                }

                operator const State_Container &() const { return cur_state; }

            public:
                State_Container cur_state;
            };
        } // Coefficients
    } // EqSolver
} // GPN
