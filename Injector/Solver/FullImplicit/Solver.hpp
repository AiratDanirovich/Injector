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

                using cRHS_t = const Eigen::VectorX<RealType>;
                using RHS_t = Eigen::VectorX<RealType>;

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

                template <typename Coefs_t>
                struct EquationView
                {
                    using MatrixRow_t = Eigen::Block<Eigen::SparseMatrix<RealType>, 1, -1, false>;
                    EquationView(MatrixRow_t A,
                                 RealType &rhs,
                                 const ptrdiff_t diag_id,
                                 const Coefs_t &neib_ids)
                        : matrix_row{A}, rhs{rhs},
                          diag_id{diag_id},
                          neib_ids{neib_ids}
                    {
                        assert(matrix_row.cols() > 1ll);
                        for (const auto id : neib_ids)
                        {
                            assert(id >= 0ll);
                            assert(id < matrix_row.cols()*matrix_row.cols());
                        }
                    }

                    MatrixRow_t matrix_row;
                    RealType &rhs;
                    const ptrdiff_t diag_id;
                    const Coefs_t &neib_ids;

                    void set_type_I(const RealType val)
                    {
                        // set diagonal value = 1.0
                        matrix_row.coeffRef(diag_id) = 1.0;
                        // set non-diagonal values = 0.0
                        for (const auto id : neib_ids)
                            matrix_row.coeffRef(id) = 0.0;
                        //  set rhs = val to satisfy: 1.0*T = val
                        rhs = val;
                    }
                    void add_rhs_type_II(const RealType val)
                    {
                        bool flag{matrix_row.coeffRef(diag_id) == 1.0};
                        for (const auto id : neib_ids)
                            flag = flag && (matrix_row.coeffRef(id) == 0.0);
                        // add the given flux to the rhs,
                        // if type_I BC was not applied 
                        // from the other face
                        if(flag == false)
                            rhs += val;
                    }
                };

                auto advance(RealType tau)
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
                    const Eigen::ArrayX<RealType> tau_factor{
                        time_factor.Divide(tau).reshaped(A_size, 1ll).eval()};

                    assert(A_size == tau_factor.size());

                    assemble_y(tripletList);
                    assemble_x(tripletList);

                    A.setFromTriplets(tripletList.begin(), tripletList.end());
                    A.diagonal() = A.diagonal() + tau_factor.matrix();

                    RHS_t rhs{
                        (state.cur_state.reshaped(A_size, 1ll).array() * tau_factor).matrix()};
                    // BC
                    applyBC(A, rhs);

                    state.cur_state = solve_linear_problem(A, rhs).array().reshaped(first_coord_size, second_coord_size);

                    cur_time += tau;

                    return std::pair{std::move(A), std::move(rhs)};
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
                    for (auto row{0ll}; row < first_coord_size; ++row)
                    {
                        // copy Laplace term in y-direction for a fixed x
                        const SpMatrix &A{splitX.LaplaceTerm(row)};
                        // convection term
                        // exclude rightmost edge
                        const auto &flow{split_flow_field.row(row).head(second_coord_size).matrix().transpose()};

                        // upper diagonal
                        for (auto col{0ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + first_coord_size, A.coeff(col, col + 1ll));
                        }

                        // main diagonal
                        const auto diag{(A.diagonal() + flow).eval()};
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
                                A.coeff(col, col - 1ll) - flow(col));
                        }
                    }
                }

                void assemble_y(auto &tripletList)
                {
                    const auto &split_flow_field{
                        convection_factory.get_flow_in_axes1()};

                    // take every line for a fixed y-node.
                    // It is a col of 2D grid representation
                    for (auto col{0ll}; col < second_coord_size; ++col)
                    {
                        // Laplace term
                        const SpMatrix &A{splitY.LaplaceTerm(col)};
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
                        for (auto row{0ll}; row < first_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            tripletList.emplace_back(l, l + 1ll, A.coeff(row, row + 1ll) + flow_minus(row + 1ll));
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
                    if (lu.info() != Eigen::Success)
                    {
                        throw std::runtime_error("LU decomposition failed!");
                    }
                    return lu.solve(b);
                }

                void applyBC(SpMatrix &A, Eigen::VectorX<RealType> &b)
                {
                    { // west face
                        const auto col{0ll};
                        {
                            auto row{0ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 2ull>{l + 1, l + first_coord_size}};
                            bc.set_west_val(view, row);
                        }
                        for (auto row{1ll}; row < first_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 3ull>{l - 1, l + 1, l + first_coord_size}};
                            bc.set_west_val(view, row);
                        }
                        {
                            auto row{first_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 2ull>{l - 1, l + first_coord_size}};
                            bc.set_west_val(view, row);
                        }
                    }
                    { // east face
                        const auto col{second_coord_size - 1ll};
                        {
                            auto row{0ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 2ull>{l + 1, l - first_coord_size}};
                            bc.set_east_val(view, row);
                        }
                        for (auto row{1ll}; row < first_coord_size - 1ll; ++row)
                        {
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 3ull>{l - 1, l + 1, l - first_coord_size}};
                            bc.set_east_val(view, row);
                        }
                        {
                            auto row{first_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 2ull>{l - 1, l - first_coord_size}};
                            bc.set_east_val(view, row);
                        }
                    }

                    { // south face
                        const auto row{0ll};
                        {
                            const auto col{0ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 2ll>{l + 1, l + first_coord_size}};
                            bc.set_south_val(view, col);
                        }
                        for (auto col{1ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 3ll>{l - first_coord_size, l + 1, l + first_coord_size}};
                            bc.set_south_val(view, col);
                        }
                        {
                            const auto col{second_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 2ll>{l - first_coord_size, l + 1}};
                            bc.set_south_val(view, col);
                        }
                    }
                    { // north face
                        const auto row{first_coord_size - 1ll};
                        {
                            const auto col{0ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 2ll>{l - 1, l + first_coord_size}};
                            bc.set_north_val(view, col);
                        }
                        for (auto col{1ll}; col < second_coord_size - 1ll; ++col)
                        {
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 3ll>{l - first_coord_size, l - 1, l + first_coord_size}};
                            bc.set_north_val(view, col);
                        }
                        {
                            const auto col{second_coord_size - 1ll};
                            const auto l{grid->to_linear(row, col)};
                            EquationView view{
                                A.row(l), b(l),
                                l,
                                std::array<ptrdiff_t, 3ll>{l - first_coord_size, l - 1}};
                            bc.set_north_val(view, col);
                        }
                    }
                }
            };
        } // SplittingMethod
    } // EqSolver
} // GPN