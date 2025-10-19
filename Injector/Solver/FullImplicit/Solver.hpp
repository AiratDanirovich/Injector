#pragma once

#include <vector>
#include <array>
#include <tuple>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/Grids1D.hpp>

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

                using cRHS_t = const Eigen::VectorX<RealType>;
                using RHS_t = Eigen::VectorX<RealType>;

                template <
                    typename LaplaceFactor_t>
                Solver(
                    const LaplaceFactor_t &laplace_factor,
                    const cptr<Grid_t> grid,
                    const Capacity_t &time_factor,
                    ptr<ConvectionTermFactory_t> convection_factory,
                    const State::State2D &initial_state,
                    const BC_t &bc,
                    const RealType initial_moment)
                    : splitX{laplace_factor, grid},
                      splitY{laplace_factor, grid},
                      time_factor{time_factor, *grid},
                      grid{grid},
                      convection_factory{convection_factory},
                      first_coord_size{grid->first_coord().mesh_size()},
                      second_coord_size{grid->second_coord().mesh_size()},
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
        //        Solver(const Solver &) = default;
        //        Solver(Solver &&) noexcept = default;

                void save_state()
                {
                    time_moments.push_back(cur_time);
                    states.emplace_back(state);
                }

                template <typename Coefs_t>
                struct EquationView
                {
                    using Matrix_t = Eigen::SparseMatrix<RealType>; // Eigen::Block<Eigen::SparseMatrix<RealType>, 1, -1, false>;
                    EquationView(Matrix_t &A,
                                 RealType &rhs,
                                 const ptrdiff_t diag_id,
                                 const Coefs_t &neib_ids)
                        : matrix{A}, rhs{rhs},
                          diag_id{diag_id},
                          neib_ids{neib_ids}
                    {
                        assert(matrix.cols() > 1ll);
                        for (const auto id : neib_ids)
                        {
                            assert(id >= 0ll);
                            assert(id < matrix.cols() * matrix.cols());
                        }
                    }

                    Matrix_t &matrix;
                    RealType &rhs;
                    const ptrdiff_t diag_id;
                    const Coefs_t &neib_ids;

                    void set_type_I(const RealType val)
                    {
                        // set diagonal value = 1.0
                        matrix.coeffRef(diag_id, diag_id) = 1.0;
                        // set non-diagonal values = 0.0
                        for (const auto id : neib_ids)
                            matrix.coeffRef(diag_id, id) = 0.0;
                        //  set rhs = val to satisfy: 1.0*T = val
                        rhs = val;
                    }
                    void add_rhs_type_II(const RealType val)
                    {
                        bool flag{matrix.coeffRef(diag_id, diag_id) == 1.0};
                        for (const auto id : neib_ids)
                            flag = flag && (matrix.coeffRef(diag_id, id) == 0.0);
                        // add the given flux to the rhs,
                        // if type_I BC was not applied
                        // from the other face
                        if (flag == false)
                            rhs += val;
                    }
                };

                auto advance(RealType tau)
                {
                    // update convection field
                    convection_factory->set_flow_field(cur_time, tau);
                    // update types of boundary conditions
                    bc.set_bc_type(cur_time+tau);

                    ptrdiff_t A_size{first_coord_size * second_coord_size};
                    assert(A_size == grid->mesh_size());
                    Eigen::SparseMatrix<RealType> A{// ctor for matrix
                                                    A_size,
                                                    A_size};

                    A.reserve(A_size * 5ll);

                    std::vector<Eigen::Triplet<RealType, ptrdiff_t>> tripletList;
                    tripletList.reserve(6ll * A_size);

                    // tau_factor = capacity/tau multiplies Delta_u at different time moments,
                    // i.e., t and t+tau
                    const Eigen::ArrayXX<RealType> tau_factor{
                        time_factor.Divide(tau).eval()};

                    assert(A_size == tau_factor.size());

                    assemble_y(tripletList);
                    assemble_x(tripletList);

                    A.setFromTriplets(tripletList.begin(), tripletList.end());
                    A.diagonal() = A.diagonal() + tau_factor.reshaped(A_size, 1ll).matrix();

                    RHS_t rhs{assemble_RHS(state, tau_factor, A_size)};
                    // BC
                    applyBC(A, rhs);

                    state.cur_state = solve_linear_problem(A, rhs).array().reshaped(first_coord_size, second_coord_size);

                    cur_time += tau;

                    return std::pair{std::move(A), std::move(rhs)};
                }

                RHS_t assemble_RHS(
                    const auto state,
                    const auto tau_factor,
                    const auto A_size) const
                {
                    return 
                        (state.cur_state.array() * tau_factor +
                        convection_factory->get_spatial_JT_contribution()).reshaped(A_size, 1ll)
                            .matrix();
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
                    const auto &split_flow_field_pos{
                        convection_factory->get_heat_flow_in_axes2_pos()};
                    const auto &split_flow_field_neg{
                        convection_factory->get_heat_flow_in_axes2_neg()};

                    const auto &pressure{convection_factory->get_pressure_field()};

                    //  Take every line for a fixed x node.
                    //  It is a row of 2D grid representation
                    for (auto row{0ll}; row < first_coord_size; ++row)
                    {
                        // copy Laplace term in y-direction for a fixed x
                        const SpMatrix &A{splitX.LaplaceTerm(row)};
                        const auto &flow_plus{split_flow_field_pos.row(row)};
                        const auto &flow_minus{split_flow_field_neg.row(row)};

                        // upper diagonal
                        for (auto col{0ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + first_coord_size,
                                                     A.coeff(col, col + 1ll) + flow_minus(col + 1ll));
                        }

                        // main diagonal
                        const auto diag{(
                                            A.diagonal() +
                                            (flow_plus.matrix().head(second_coord_size) -
                                             flow_minus.matrix().tail(second_coord_size))
                                                .transpose())
                                            .eval()};

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
                                A.coeff(col, col - 1ll) - flow_plus(col));
                        }
                    }
                }

                void assemble_y(auto &tripletList)
                {
                    const auto &split_flow_field_pos{
                        convection_factory->get_heat_flow_in_axes1_pos()};
                    const auto &split_flow_field_neg{
                        convection_factory->get_heat_flow_in_axes1_neg()};

                    // take every line for a fixed y-node.
                    // It is a col of 2D grid representation
                    for (auto col{0ll}; col < second_coord_size; ++col)
                    {
                        // Laplace term
                        const SpMatrix &A{splitY.LaplaceTerm(col)};

                        const auto &flow_plus{split_flow_field_pos.col(col)};
                        const auto &flow_minus{split_flow_field_neg.col(col)};

                        // upper diagonal
                        for (auto row{0ll}; row < first_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + 1ll,
                                                     A.coeff(row, row + 1ll) + flow_minus(row + 1ll));
                        }

                        // main diagonal
                        const auto diag{(
                                            A.diagonal() +
                                            flow_plus.matrix().head(first_coord_size) -
                                            flow_minus.matrix().tail(first_coord_size))
                                            .eval()};
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
                                A.coeff(row, row - 1ll) - flow_plus(row));
                        }
                    }
                }

                const SplittingMethod::SplitX splitX;
                const SplittingMethod::SplitY splitY;
                BC_t bc;
                const TemporalTerm time_factor;
                double cur_time;
                // by reference!
                ptr<ConvectionTermFactory_t> convection_factory;
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
        } // SplittingMethod
    } // EqSolver
} // GPN