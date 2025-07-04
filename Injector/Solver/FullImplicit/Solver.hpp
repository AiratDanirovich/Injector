#pragma once

#include <vector>
#include <tuple>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids1D.hpp>
#include <Injector/Properties/PhysicalField.hpp>

#include <Injector/Solver/CapacityTerm.hpp>
#include <Injector/Solver/SplittingMethod/SplitX.hpp>
#include <Injector/Solver/SplittingMethod/SplitY.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace FullImplicit
        {
            template <
                typename Grid_t,
                typename Capacity_t,
                typename ConvectionTermFactory_t,
                typename BC_t>
            struct Solver
            {
                using SpMatrix = SplittingMethod::SpMatrix;

                using RHS_t = Eigen::Map<const Eigen::VectorX<RealType>>;

                template <
                    typename LaplaceFactor_t>
                Solver(
                    const LaplaceFactor_t &laplace_factor,
                    const cptr<Grid_t> grid,
                    const Capacity_t &time_factor,
                    ConvectionTermFactory_t &convection_factory,
                    const State::State2D &initial_state,
                    const BC_t &bc,
                    RealType initial_moment)
                    : splitX{laplace_factor, grid},
                      splitY{laplace_factor, grid},
                      time_factor{time_factor, *grid},
                      grid{grid},
                      convection_factory{convection_factory},
                      first_coord_size{grid->first_coord.mesh_size()},
                      second_coord_size{grid->second_coord.mesh_size()},
                      state{initial_state}, // init with initial condition
                      bc{bc},
                      states(),
                      time_moments(),
                      cur_time{initial_moment}
                {
                    time_moments.reserve(10ull);
                    states.reserve(10ull);
                    save_state();
                }
                Solver(const Solver &) = default;
                Solver(Solver &&) noexcept = default;

                void save_state()
                {
                    time_moments.push_back(cur_time);
                    states.emplace_back(state);
                }

                void advance(RealType tau)
                {
                    // update convection field
                    convection_factory.set_flow_field(cur_time, tau);
                    // update boundary conditions
                    bc.set_vals(cur_time + tau);

                    ptrdiff_t A_size{first_coord_size * second_coord_size};
                    assert(A_size == grid->mesh_size());
                    Eigen::SparseMatrix<RealType> A{// ctor per matrix
                                                    A_size,
                                                    A_size};
                    A.reserve(A_size * 5ll);

                    std::vector<Eigen::Triplet<RealType, ptrdiff_t>> tripletList;
                    tripletList.reserve(5ll * A_size);

                    // tau_factor multiplies Delta_u at different time moments,
                    // i.e., t and t+tau
                    const auto tau_factor{
                        time_factor.Divide(tau).eval()};

                    assert(A_size == tau_factor.size());

                    //        assemble_y(tripletList);
                    //        assemble_x(tripletList);

                    A.setFromTriplets(tripletList.begin(), tripletList.end());
                    const RHS_t capacity_term(tau_factor.data(), A_size, 1ll);
                    A.diagonal() = A.diagonal() + capacity_term;

                    Eigen::Map<Eigen::MatrixXf> rhs{
                        (state.cur_state.array() * time_factor.array()).matrix()};

                    // BC
                    applyBC_split_x(A, rhs, 0ll);
                    applyBC_split_y(A, rhs, 0ll);

                    const auto val{solve_linear_problem(A, rhs)};

                    cur_time += tau;
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
                void assemble_x(auto &tripletList)
                {
                    const auto &split_flow_field{
                        convection_factory.get_flow_in_axes2()};

                    //  Take every line for a fixed x node.
                    //  It is a row of 2D grid representation
                    for (std::ptrdiff_t row = 0; row < first_coord_size; ++row)
                    {
                        // copy Laplace term in y-direction for a fixed x
                        const SpMatrix &A{splitX.LaplaceTerm(row)};
                        // convection term
                        // exclude rightmost edge
                        const auto &flow{split_flow_field.row(row).head(second_coord_size).matrix().transpose()};

                        // upper diagonal
                        for (std::ptrdiff_t col{1ll}; col < second_coord_size; ++col)
                        {
                            const auto node_id{};
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + first_coord_size, A.coeff(col - 1ll, col));
                        }

                        // main diagonal
                        const auto diag{(A.diagonal() + flow).eval()};
                        for (std::ptrdiff_t col{0ll}; col < second_coord_size; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l, diag(col));
                        }

                        // lower diagonal
                        for (std::ptrdiff_t col{0ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(
                                l, l - first_coord_size,
                                A.coeff(col + 1ll, col) - flow(col + 1ll));
                        }
                    }
                }

                void assemble_y(auto &tripletList)
                {
                    const auto &split_flow_field{
                        convection_factory.get_flow_in_axes1()};

                    // take every line for a fixed y-node.
                    // It is a col of 2D grid representation
                    for (std::ptrdiff_t col = 0; col < second_coord_size; ++col)
                    {
                        // Laplace term
                        const SpMatrix& A{splitY.LaplaceTerm(col)};
                        // convection term
                        const auto &temp_flow{split_flow_field.col(col).matrix()};
                        // positive flow values
                        const auto flow_plus{(temp_flow.array() + temp_flow.array().abs()) / 2.0};
                        assert(flow_plus.rows() == first_coord_size + 1ll);
                        assert(std::all_of(flow_plus.cbegin(), flow_plus.cend(), [](const RealType v)
                                           { return v >= 0.0; }));
                        // negative flow values
                        const auto flow_minus{(temp_flow.array() - temp_flow.array().abs()) / 2.0};
                        assert(flow_minus.rows() == first_coord_size + 1ll);
                        assert(std::all_of(flow_minus.cbegin(), flow_minus.cend(), [](const RealType v)
                                           { return v <= 0.0; }));

                        // upper diagonal
                        for (std::ptrdiff_t row{1ll}; row < second_coord_size; ++row)
                        {
                            const auto node_id{};
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + 1ll, A.coeff(row - 1ll, row) + flow_minus(row) );
                        }

                        // main diagonal
                        const auto diag{(
                            A.diagonal() + 
                            flow_plus.matrix().head(first_coord_size) - 
                            flow_minus.matrix().tail(first_coord_size)).eval()};
                        for (std::ptrdiff_t row{0ll}; row < second_coord_size; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l, diag(row));
                        }

                        // lower diagonal
                        for (std::ptrdiff_t row{0ll}; row < second_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(
                                l, l - 1ll,
                                A.coeff(row + 1ll, row) - flow_plus(row + 1ll));
                        }
                    }
                }

                const SplittingMethod::SplitX splitX;
                const SplittingMethod::SplitY splitY;
                BC_t bc;
                const TemporalTerm time_factor;
                double cur_time;
                // by reference!
                ConvectionTermFactory_t &convection_factory;
                // required to keep grid in memory ////
                const cptr<Grid_t> grid; //////////////
                ///////////////////////////////////////
                const std::ptrdiff_t first_coord_size;
                const std::ptrdiff_t second_coord_size;
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
                    {
                        BoundaryConditions::MatrixView view{A.row(0ll), b.row(0ll), 0ll, 1ll};
                        bc.set_west_val(view, i);
                    }
                    {
                        std::ptrdiff_t n = A.outerSize() - 1;
                        BoundaryConditions::MatrixView view{A.row(n), b.row(n), n, n - 1};
                        bc.set_east_val(view, i);
                    }
                }

                void applyBC_split_y(SpMatrix &A, RHS_t &b, ptrdiff_t j)
                {
                    {
                        BoundaryConditions::MatrixView view{A.row(0ll), b.row(0ll), 0ll, 1ll};
                        bc.set_south_val(view, j);
                    }
                    {
                        std::ptrdiff_t n = A.outerSize() - 1;
                        BoundaryConditions::MatrixView view{A.row(n), b.row(n), n, n - 1};
                        bc.set_north_val(view, j);
                    }
                }
            };
        } // SplittingMethod
    } // EqSolver
} // GPN