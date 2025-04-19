#pragma once

#include <vector>
#include <tuple>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/PhysicalField.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>

#include <Injector/Solver/SplittingMethod/SplitX.hpp>
#include <Injector/Solver/SplittingMethod/SplitY.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace SplittingMethod
        {
            template <typename Capacity_t>
            struct TemporalTerm
            {
                template <typename Grid_t>
                TemporalTerm(
                    const Capacity_t& factor,
                    const Grid_t &grid)
                    : factor{grid.volumes()*factor.values()} // volumes are taken into account
                {
                }

                auto Divide(RealType tau) const
                {
                    return (static_cast<RealType>(1.0) / tau) * factor;
                }

            protected:
                // multiplied by cell volume
                Eigen::ArrayXX<RealType> factor;
            };

            template <
                typename Grid_t,
                typename Capacity_t,
                typename LaplaceFactor_t>
            struct Solver
            {
                using Map1D =
                    Eigen::Map<
                        Eigen::ArrayX<RealType>>;

                using Map1D_Stride =
                    Eigen::Map<
                        Eigen::VectorX<RealType>,
                        0,
                        Eigen::OuterStride<Eigen::Dynamic>>;

                using RHS_t = Eigen::VectorX<RealType>;

                Solver(
                    const Grid_t &grid,
                    const Capacity_t &time_factor,
                    const LaplaceFactor_t &laplace_factor,
                    const State2D_t &initial_state,
                    const BC_t &bc,
                    Realtype initial_moment = 0.0)
                    : splitX{properties, grid},
                      splitY{properties, grid},
                      time_factor{time_factor},
                      grid{grid},
                      state{initial_state}, // init with initial condition
                      bc{bc},
                      states(),
                      time_moments()
                {
                    time_moments.reserve(10ull);
                    states.reserve(10ull);
                    time_moments.push_back(initial_moment);
                    states.emplace_back(state);
                }

                template <typename Factory_t>
                static auto set_from_factory(const Factory_t &factory)
                {
                    return Solver{
                        factory.grid, factory.capacity,
                        factory.laplace_factor,
                        factory.initial_state,
                        factory.bc, factory.initial_moment};
                }

                void advance(RealType tau)
                {
                    // tau_factor multiplies Delta_u at different time moments,
                    // i.e., t and t+tau
                    Eigen::ArrayXX<RealType> tau_factor{
                        factor.Divide(tau)};

                    // solve a set of 1D problems in y-direction, for various x-coords
                    bc.set_vals(time_moments.back() + tau/2.0);

                    solve_split_x(tau_factor.data());

                    // solve a set of 1D problems in x-direction, for various y-coords
                    bc.set_vals(time_moments.back() + tau);
                    solve_split_y(tau_factor.data());

                    time_moments.push_back(time_moments.back() + tau);
                    states.emplace_back(state);
                }

                struct Solution
                {
                    Solution(
                        const std::vector<RealType> &times,
                        const std::vector<State::State2D> &states)
                        : times{times}, states{states}
                    {
                    }

                    auto back() const
                    {
                        return std::pair{times.back(), states.back()};
                    }

                    const std::vector<RealType> &times;
                    const std::vector<State::State2D> &states;
                };

                Solution solution() const
                {
                    return {time_moments, states};
                }

            protected:
                void solve_split_x(RealType *tau_factor)
                {
                    // nmbr of nodes in the 1D problem
                    ptrdiff_t chunk_size = grid.second_coord.size();
                    // stride in a 1D layout of 2D unknown temperature values
                    ptrdiff_t stride_size = grid.first_coord.size();
                    // to be provided to Eigen::Map
                    Eigen::OuterStride<Eigen::Dynamic> stride{stride_size};

#pragma omp parallel for num_threads(16) schedule(dynamic)
                    //  take every line along x-direction. A line per y-node
                    for (std::ptrdiff_t i = 0; i < grid.first_coord.size(); ++i)
                    {
                        // memory chunk in capacity-container, corresponding to x-line
                        const auto time_factor =
                            Map1D_Stride{
                                tau_factor + i,
                                chunk_size, stride};

                        // memory chunk in temperature, corresponding to x-line
                        auto data =
                            Map1D_Stride{
                                state.cur_state.data() + i,
                                chunk_size, stride};

                        // right handside of Au = b problem
                        // source is only assumed in the last split step
                        RHS_t rhs =
                            (data.array() * time_factor.array()).matrix();

                        SpMatrix A{splitX.LaplaceTerm(i)};
                        A.diagonal() = A.diagonal() + time_factor;
                        applyBC_split_x(A, rhs, i);
                        // update current state
                        data = solve_linear_problem(A, rhs);
                    }
                }

                void solve_split_y(RealType *tau_factor)
                {
#pragma omp parallel for
                    // take every line along x-direction. A line per y-node.
                    // It is a col of 2D grid representation
                    for (std::ptrdiff_t j = 0; j < grid.second_coord.size(); ++j)
                    {
                        // memory chunk in capacity-container, corresponding to x-line
                        const auto time_factor =
                            Map1D{
                                tau_factor + grid.first_coord.size() * j,
                                grid.first_coord.size()};

                        // memory chunk in temperature, corresponding to x-line
                        auto data =
                            Map1D{
                                state.cur_state.data() + grid.first_coord.size() * j,
                                grid.first_coord.size()};

                        const auto source = 0.0;
                            // Map1D{
                            //     properties->source_vol_f.data() + grid.first_coord.size() * j,
                            //     (ptrdiff_t)grid.first_coord.size()};

                        // right handside of Au = b problem
                        RHS_t rhs = (data * time_factor + source).matrix();

                        SpMatrix A{splitY.LaplaceTerm(j)};
                        A.diagonal() = A.diagonal() + time_factor.matrix();
                        applyBC_split_y(A, rhs, j);
                        // update current state
                        data = solve_linear_problem(A, rhs);
                    }
                }

                SplitX splitX;
                SplitY splitY;
                BoundaryConditions::BoundaryConditions bc;
                TemporalTerm<Capacity_t> time_factor;
                const Grid_t &grid;
                State::State2D state;
                std::vector<State::State2D> states;
                std::vector<RealType> time_moments;

            protected:
                Eigen::SparseLU<SpMatrix> lu;

                Eigen::VectorXd solve_linear_problem(
                    const SpMatrix &A,
                    const Eigen::VectorX<RealType> &b)
                {
                    // A.makeCompressed();
                    Eigen::SparseLU<SpMatrix> lu;
                    lu.analyzePattern(A); // this is common for every matrix A. Can be optimized
                    lu.factorize(A);
                    return lu.solve(b);
                }

                void applyBC_split_x(SpMatrix &A, RHS_t &b, ptrdiff_t i)
                {
                    A.coeffRef(0, 0) = 1;
                    A.coeffRef(0, 1) = 0;
                    b(0) = bc.west_vals(i);

                    ptrdiff_t n = A.outerSize() - 1;
                    A.coeffRef(n, n) = 1.0;
                    A.coeffRef(n, n - 1) = 0.0;
                    b(n) = bc.east_vals(i);
                }

                void applyBC_split_y(SpMatrix &A, RHS_t &b, ptrdiff_t j)
                {
                    A.coeffRef(0, 0) = 1;
                    A.coeffRef(0, 1) = 0;
                    b(0) = bc.south_vals(j);

                    ptrdiff_t n = A.outerSize() - 1;
                    A.coeffRef(n, n) = 1.0;
                    A.coeffRef(n, n - 1) = 0.0;
                    b(n) = bc.north_vals(j);
                }
            };
        } // SplittingMethod
    } // EqSolver
} // GPN