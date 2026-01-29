/**
 * @file AppendSolver.hpp
 * @brief Solves the reservoir flow problem with the rate control at bottomhole
 *
 * Originally it is a copy-paste of Solver/FullImplicit/Solver.hpp.
 * In the nearest future they must be combined into a unified solver-library
 * with specialization between various types of problems
 *
 * @author Arthur Salamatin
 * @date 2026-01-28
 */

#pragma once

#include <vector>
#include <array>
#include <tuple>
#include <cassert>
#include <type_traits>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids1D.hpp>

#include <Injector/Solver/CapacityTerm.hpp>
#include <Injector/Solver/SplittingMethod/SplitX.hpp>
#include <Injector/Solver/SplittingMethod/SplitY.hpp>
#include <Injector/Solver/EquationView.hpp>

#include <Injector/Model/Well/WellBottomHoleRateControl.hpp>

namespace GPN
{
    namespace EqSolver
    {
        using TripletContainer = std::vector<Eigen::Triplet<RealType, ptrdiff_t>>;

        namespace FullImplicit
        {
            template <
                typename Grid2D_t,
                typename Capacity_t,
                typename ConvectionTermFactory_t,
                typename BC_t,
                typename Well_t>
            struct WellRateControlSolver
            {
                using SpMatrix = SplittingMethod::SpMatrix;

                using cRHS_t = const Eigen::VectorX<RealType>;
                using RHS_t = Eigen::VectorX<RealType>;

                template <
                    typename LaplaceFactor_t>
                WellRateControlSolver(
                    const LaplaceFactor_t &laplace_factor,
                    const cptr<Grid2D_t> grid,
                    const Capacity_t &time_factor,
                    ptr<ConvectionTermFactory_t> convection_factory,
                    const State::State2D &initial_state,
                    const BC_t &bc,
                    const Well_t &well,
                    const RealType initial_moment,
                    const RealType P_bot_init = 0.0)
                    : splitX{laplace_factor, grid},
                      splitY{laplace_factor, grid},
                      time_factor{time_factor, *grid},
                      grid{grid},
                      convection_factory{convection_factory},
                      first_coord_size{grid->first_coord().mesh_size()},
                      second_coord_size{grid->second_coord().mesh_size()},
                      state{initial_state}, // init with initial condition
                      bc{bc},
                      well{well},
                      states{},
                      time_moments{},
                      P_bots{},
                      P_bot_memory{P_bot_init},
                      cur_time{initial_moment}
                {
                    time_moments.reserve(10ull);
                    states.reserve(10ull);
                    save_state();
                }

                void save_state()
                {
                    time_moments.push_back(cur_time);
                    states.emplace_back(state);
                    P_bots.emplace_back(P_bot);
                }

                TripletContainer set_triplets(const ptrdiff_t A_size, const RealType tau)
                {
                    TripletContainer tripletList;

                    assert(A_size == grid->mesh_size() + 1ll);
                    tripletList.reserve(8ll * A_size);
#pragma region MASS-TRANSFER-IN-RESERVOIR
#pragma region TEMPORAL-CONTRIBUTION-IN-RESERVOIR
                    // take reservoir capacity into account
                    const Eigen::ArrayXX<RealType> tau_factor{
                        time_factor.Divide(tau).eval()};
                    assert(tau_factor.rows() == A_size - 1ll);
                    assert(tau_factor.cols() == A_size - 1ll);
                    // append capacity to the matri diagonal
                    for (auto row{0ll}; row < tau_factor.rows(); ++row)
                        for (auto col{0ll}; col < tau_factor.cols(); ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l, tau_factor(row, col));
                        }
#pragma endregion
#pragma region LAPLACE-TERM-IN-RESERVOIR
                    // take Laplace contribution into account
                    assemble_y_noconvection(tripletList);
                    assemble_x_noconvection(tripletList);
#pragma endregion
#pragma endregion
#pragma region WELL-CONDITION
                    // take well condition into account
                    const auto &PI{well.PI}; // well productivity index
                    // set the P_bot coefficient
                    const auto matrix_row{A_size - 1ll}; // id of unknown bottomhole pressure
                    // set last matrix diag element
                    tripletList.emplace_back(matrix_row, matrix_row, -PI.sum());
                    // set P_i,3 coeffs if non-zero
                    for (auto row{0}; row < PI.rows(); ++row)
                    { // loop over every row of the sandface
                        if (PI(row) != 0.0)
                        {
                            // set last matrix row
                            const auto matrix_col{grid->to_linear(row, 0ll)};
                            tripletList.emplace_back(matrix_row, matrix_col, PI(row));
                            // set last matrix col
                            tripletList.emplace_back(matrix_col, matrix_row, -PI(row));
                        }
                    }
#pragma endregion
                    return tripletList;
                }

                auto advance(const RealType tau)
                {
                    ptrdiff_t A_size{first_coord_size * second_coord_size + 1ll};
                    Eigen::SparseMatrix<RealType> A{// ctor for matrix, reservoir + bottomwell pressure
                                                    A_size,
                                                    A_size};
                    A.reserve(A_size * 6ll);

                    TripletContainer tripletList{get_triplets(A_size, tau)};
                    assert(std::all_of(
                        tripletList.cbegin(), 
                        tripletList.cend(), 
                        [](const auto& v){
                            return !std::isnan(v.value()) && !std::isinf(v.value());
                        }));
                    A.setFromTriplets(tripletList.begin(), tripletList.end());

                    const Eigen::ArrayXX<RealType> tau_factor{
                        time_factor.Divide(tau).eval()};
                    RHS_t rhs{assemble_RHS_noconvection(state, tau_factor, A_size)};
                    
                    assert(std::all_of(
                        rhs.cbegin(), 
                        rhs.cend(), 
                        [](const auto& v){
                            return !std::isnan(v) && !std::isinf(v);
                        }));

                    // BC
                    // update types of boundary conditions
                    bc.set_bc_type(cur_time + tau);
                    applyBC(A, rhs);

                    const auto solution{solve_linear_problem(A, rhs).array().reshaped(first_coord_size, second_coord_size)};
                    
                    for(auto it{solution.cbegin()}; it != solution.cend(); ++it)
                        {
                            assert(!std::isnan(*it) && !std::isinf(*it));
                        }

                    state.cur_state = solution.topRows(A_size-2ll);
                    P_bot_memory = solution.bottomRows(1ll);

                    cur_time += tau;

                    return std::pair{std::move(A), std::move(rhs)};
                }

                RHS_t assemble_RHS_noconvection(
                    const auto state,
                    const auto tau_factor,
                    const auto A_size) const
                {
                    RHS_t rhs{RHS_t::Zero(A_size)};
                    rhs.topRows(A_size-1ll) = (state.cur_state.array() * tau_factor)
                                          .reshaped(A_size, 1ll)
                                          .matrix();
                    rhs.bottomRows(1ll) = well.history->rate();
                    return rhs;
                }

                struct Solution
                {
                    Solution(
                        const std::vector<RealType> &times,
                        const std::vector<State::State2D> &states,
                        const std::vector<RealType> & P_bot)
                        : times{times}, states{states}, P_bot{P_bot}
                    {
                    }

                    auto back() const
                    {
                        return std::pair{times.back(), states.back()};
                    }

                    const std::vector<RealType> &times;
                    const std::vector<RealType> &P_bot;
                    const std::vector<State::State2D> &states;
                };

                Solution solution() const
                {
                    return {time_moments, states, P_bot};
                }

                const State::State2D &get_state() const
                {
                    return state;
                }

                const RealType P_bot() const
                {
                    return P_bot_memory;
                }

            protected:
                void assemble_x_noconvection(auto &tripletList)
                {
                    //  Take every line for a fixed x node.
                    //  It is a row of 2D grid representation
                    for (auto row{0ll}; row < first_coord_size; ++row)
                    {
                        // copy Laplace term in y-direction for a fixed x
                        const SpMatrix &A{splitX.LaplaceTerm(row)};

                        // upper diagonal
                        for (auto col{0ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + first_coord_size,
                                                     A.coeff(col, col + 1ll));
                        }

                        // main diagonal
                        const auto diag{(A.diagonal()).eval()};
                        assert(diag.size() == second_coord_size);

                        for (auto col{0ll}; col < second_coord_size; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l, diag(col));
                        }

                        // lower diagonal
                        for (auto col{1ll}; col < second_coord_size; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            assert(l >= first_coord_size);
                            tripletList.emplace_back(
                                l, l - first_coord_size,
                                A.coeff(col, col - 1ll));
                        }
                    }
                }

                void assemble_y_noconvection(auto &tripletList)
                {
                    // take every line for a fixed y-node.
                    // It is a col of 2D grid representation
                    for (auto col{0ll}; col < second_coord_size; ++col)
                    {
                        // Laplace term
                        const SpMatrix &A{splitY.LaplaceTerm(col)};

                        // upper diagonal
                        for (auto row{0ll}; row < first_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + 1ll,
                                                     A.coeff(row, row + 1ll));
                        }

                        // main diagonal
                        const auto diag{(A.diagonal()).eval()};
                        assert(diag.size() == first_coord_size);
                        for (auto row{0ll}; row < first_coord_size; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l, diag(row));
                        }

                        // lower diagonal
                        for (auto row{1ll}; row < first_coord_size; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            assert(l >= 1ll);
                            tripletList.emplace_back(
                                l, l - 1ll,
                                A.coeff(row, row - 1ll));
                        }
                    }
                }

                const SplittingMethod::SplitX splitX;
                const SplittingMethod::SplitY splitY;
                BC_t bc;
                const Well_t &well;
                const TemporalTerm time_factor;
                double cur_time;
                // by reference!
                ptr<ConvectionTermFactory_t> convection_factory;
                // required to keep grid in memory ////
                const cptr<Grid2D_t> grid; //////////////
                ///////////////////////////////////////
                const std::ptrdiff_t first_coord_size;
                const std::ptrdiff_t second_coord_size;
                State::State2D state;
                RealType P_bot_memory;
                std::vector<State::State2D> states;
                std::vector<RealType> time_moments, P_bots;

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
                    if (lu.info() != Eigen::Success)
                    {
                        throw std::runtime_error("LU decomposition failed!");
                    }
                    return lu.solve(b);
                }

                /// @brief Apply boundary conditions of the problem
                /// @param A NxN sparse matrix of the problem with a row per every node of the 2D mesh
                /// @param b Nx1 vector of the rhs of the algebraic problem
                void applyBC(SpMatrix &A, Eigen::VectorX<RealType> &b)
                {
                    { // west face
                        const auto col{0ll};
                        {
                            auto row{0ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neibs{std::array<ptrdiff_t, 2ull>{l + 1, l + first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neibs};
                            bc.set_west_val(view, row);
                        }
                        for (auto row{1ll}; row < first_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            const auto neibs{std::array<ptrdiff_t, 3ull>{l - 1, l + 1, l + first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neibs};
                            bc.set_west_val(view, row);
                        }
                        {
                            auto row{first_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neibs{std::array<ptrdiff_t, 2ull>{l - 1, l + first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neibs};
                            bc.set_west_val(view, row);
                        }
                    }
                    { // east face
                        const auto col{second_coord_size - 1ll};
                        {
                            auto row{0ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neibs{std::array<ptrdiff_t, 2ull>{l + 1, l - first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neibs};
                            bc.set_east_val(view, row);
                        }
                        for (auto row{1ll}; row < first_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            const auto neibs{std::array<ptrdiff_t, 3ull>{l - 1, l + 1, l - first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neibs};
                            bc.set_east_val(view, row);
                        }
                        {
                            auto row{first_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neibs{std::array<ptrdiff_t, 2ull>{l - 1, l - first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neibs};
                            bc.set_east_val(view, row);
                        }
                    }

                    { // south face
                        const auto row{0ll};
                        {
                            const auto col{0ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neib{std::array<ptrdiff_t, 2ll>{l + 1, l + first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neib};
                            bc.set_south_val(view, col);
                        }
                        for (auto col{1ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            const auto neib{std::array<ptrdiff_t, 3ll>{l - first_coord_size, l + 1, l + first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neib};
                            bc.set_south_val(view, col);
                        }
                        {
                            const auto col{second_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neib{std::array<ptrdiff_t, 2ll>{l - first_coord_size, l + 1}};

                            EquationView view{
                                A, b(l),
                                l,
                                neib};
                            bc.set_south_val(view, col);
                        }
                    }
                    { // north face
                        const auto row{first_coord_size - 1ll};
                        {
                            const auto col{0ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neib{std::array<ptrdiff_t, 2ll>{l - 1, l + first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neib};
                            bc.set_north_val(view, col);
                        }
                        for (auto col{1ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            const auto neib{std::array<ptrdiff_t, 3ll>{l - first_coord_size, l - 1, l + first_coord_size}};
                            EquationView view{
                                A, b(l),
                                l,
                                neib};
                            bc.set_north_val(view, col);
                        }
                        {
                            const auto col{second_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            const auto neib{std::array<ptrdiff_t, 3ll>{l - first_coord_size, l - 1}};
                            EquationView view{
                                A, b(l),
                                l,
                                neib};
                            bc.set_north_val(view, col);
                        }
                    }
                }
            };
        } // FullImplicit
    } // EqSolver
} // GPN