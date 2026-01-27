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
#include <Injector/Solver/EquationView.hpp>

namespace GPN
{
    namespace EqSolver
    {
        namespace SplittingMethod
        {
            template <
                typename Grid_t,
                typename Capacity_t,
                typename ConvectionTermFactory_t>
            struct Solver
            {
                using BC_t = BoundaryConditions::GeneralBC;


                using Map1D =
                    Eigen::Map<
                        Eigen::ArrayX<RealType>>;
                using cMap1D =
                    Eigen::Map<
                        const Eigen::ArrayX<RealType>>;

                using Stride_t = Eigen::InnerStride<Eigen::Dynamic>;

                using cMap1D_Stride =
                    Eigen::Map<
                        const Eigen::VectorX<RealType>,
                        Eigen::Unaligned,
                        Stride_t>;

                using Map1D_Stride =
                    Eigen::Map<
                        Eigen::VectorX<RealType>,
                        Eigen::Unaligned,
                        Stride_t>;

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

                void save_state()
                {
                    time_moments.push_back(cur_time);
                    states.emplace_back(state);
                }

                Solver(const Solver &) = default;
                Solver(Solver &&) noexcept = default;

                template <typename Factory_t>
                static auto set_from_factory(const Factory_t &factory)
                {
                    return Solver{
                        factory.grid, factory.capacity,
                        factory.laplace_factor,
                        factory.initial_state,
                        factory.bc, factory.initial_moment};
                }

                // template <typename Coefs_t>
                // struct EquationView
                // {
                //     using Matrix_t = Eigen::SparseMatrix<RealType>;
                //     EquationView(Matrix_t &A,
                //                  RealType &rhs,
                //                  const ptrdiff_t diag_id,
                //                  const Coefs_t &neib_ids)
                //         : matrix{A}, rhs{rhs},
                //           diag_id{diag_id},
                //           neib_ids{neib_ids}
                //     {
                //         assert(matrix.cols() > 1ll);
                //         for (const auto id : neib_ids)
                //         {
                //             assert(id >= 0ll);
                //             assert(id < matrix.cols() * matrix.cols());
                //         }
                //     }

                //     Matrix_t &matrix;
                //     RealType &rhs;
                //     const ptrdiff_t diag_id;
                //     const Coefs_t &neib_ids;

                //     void set_type_I(const RealType val)
                //     {
                //         // set diagonal value = 1.0
                //         matrix.coeffRef(diag_id, diag_id) = 1.0;
                //         // set non-diagonal values = 0.0
                //         for (const auto id : neib_ids)
                //             matrix.coeffRef(diag_id, id) = 0.0;
                //         //  set rhs = val to satisfy: 1.0*T = val
                //         rhs = val;
                //     }
                //     void add_rhs_type_II(const RealType val)
                //     {
                //         // add the given flux to the rhs
                //         rhs += val;
                //     }
                // };

                void advance(RealType tau)
                {
                    convection_factory.set_flow_field(cur_time, tau);
                    //    convection_factory.set_boundary_conditions(t0, t1);

                    // tau_factor multiplies Delta_u at different time moments,
                    // i.e., t and t+tau
                    Eigen::ArrayXX<RealType> tau_factor{
                        time_factor.Divide(tau / 2.0)};

// assert(false && "Fix set_vals method");


                    // solve a set of 1D problems in y-direction, for various x-coords
                //    bc.set_vals(time_moments.back() + tau / 4.0);
                    solve_split_y(tau_factor.data());

                    // solve a set of 1D problems in x-direction, for various y-coords
                //    bc.set_vals(time_moments.back() + tau / 2.0);
                    solve_split_x(tau_factor.data());

                    // solve a set of 1D problems in x-direction, for various y-coords
                //    bc.set_vals(time_moments.back() + tau * 3.0 / 4.0);
                    solve_split_x(tau_factor.data());

                    // solve a set of 1D problems in y-direction, for various x-coords
                //    bc.set_vals(time_moments.back() + tau);
                    solve_split_y(tau_factor.data());

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
                void solve_split_x(const RealType *const tau_factor)
                {
                    // nmbr of nodes in the 1D problem
                    ptrdiff_t chunk_size = second_coord_size;
                    // stride in a 1D layout of 2D unknown temperature values
                    ptrdiff_t stride_size = first_coord_size;
                    // to be provided to Eigen::Map
                    Stride_t stride{stride_size};

                    const auto &split_flow_field{
                        convection_factory.get_heat_flow_in_axes2()};

 #pragma omp parallel for // num_threads(16) schedule(dynamic)
                    //   take every line along x-direction. A line per y-node
                    //  It is a row of 2D grid representation
                    for (std::ptrdiff_t i = 0; i < first_coord_size; ++i)
                    {
                        // memory chunk in capacity-container, corresponding to x-line
                        const auto time_factor =
                            cMap1D_Stride{
                                tau_factor + i,
                                chunk_size, stride};

                        // memory chunk in temperature, corresponding to x-line
                        auto data =
                            Map1D_Stride{
                                state.cur_state.data() + i,
                                chunk_size, stride};

                        // right handside of Au = b problem
                        // source is only assumed in the last split step
                        RHS_t rhs{
                            (data.array() * time_factor.array()).matrix()};

                        // Laplace term
                        SpMatrix A{splitX.LaplaceTerm(i)};
                        // cumulative term
                        A.diagonal() = A.diagonal() + time_factor;
                        // convection term
                        const auto &temp_flow{split_flow_field.row(i).matrix()};
                        // negative flow values
                        const auto flow_plus{((temp_flow.array() + temp_flow.array().abs()) / 2.0)};
                        assert(flow_plus.cols() == second_coord_size + 1ll);
                        assert(std::all_of(flow_plus.cbegin(), flow_plus.cend(), [](const RealType v)
                                           { return v >= 0.0; }));
                        // positive flow values
                        const auto flow_minus{(temp_flow.array() - temp_flow.array().abs()) / 2.0};
                        assert(flow_minus.cols() == second_coord_size + 1ll);
                        assert(std::all_of(flow_minus.cbegin(), flow_minus.cend(), [](const RealType v)
                                           { return v <= 0.0; }));

                        // const auto &flow{split_flow_field.row(i).head(second_coord_size).matrix().transpose()};
                        // exclude rightmost edge

                        A.diagonal() = 
                            (A.diagonal() + 
                            (flow_plus.matrix().head(second_coord_size) - 
                            flow_minus.matrix().tail(second_coord_size)).transpose())
                            .eval();

                        // exclude leftmost and rightmost edges
                        for (auto idx{1ll}; idx < A.rows(); ++idx)
                            A.coeffRef(idx, idx - 1ll) -= flow_plus(idx);
                        // exclude leftmost and rightmost edges
                        for (auto idx{0ll}; idx < A.rows() - 1ll; ++idx)
                            A.coeffRef(idx, idx + 1ll) += flow_minus(idx + 1ll);

                        // BC
                        applyBC_split_x(A, rhs, i);
                        // update current state
                        data = solve_linear_problem(A, rhs);
                    }
                }

                void solve_split_y(const RealType *const tau_factor)
                {
                    const auto &split_flow_field{
                        convection_factory.get_heat_flow_in_axes1()};

#pragma omp parallel for
                    // take every line along x-direction. A line per y-node.
                    // It is a col of 2D grid representation
                    for (std::ptrdiff_t j = 0; j < second_coord_size; ++j)
                    {
                        // memory chunk in capacity-container, corresponding to x-line
                        const auto time_factor =
                            cMap1D{
                                tau_factor + first_coord_size * j,
                                first_coord_size};

                        // memory chunk in temperature, corresponding to x-line
                        auto data =
                            Map1D{
                                state.cur_state.data() + first_coord_size * j,
                                first_coord_size};

                        const auto source = 0.0;
                        // cMap1D{
                        //     properties->source_vol_f.data() + first_coord_size * j,
                        //     first_coord_size};

                        // right handside of Au = b problem
                        RHS_t rhs{(data * time_factor + source).matrix()};

                        // Laplace term
                        SpMatrix A{splitY.LaplaceTerm(j)};
                        // cumulative term
                        A.diagonal() = A.diagonal() + time_factor.matrix();
                        // convection term
                        const auto &temp_flow{split_flow_field.col(j).matrix()};
                        // negative flow values
                        const auto flow_plus{(temp_flow.array() + temp_flow.array().abs()) / 2.0};
                        assert(flow_plus.rows() == first_coord_size + 1ll);
                        assert(std::all_of(flow_plus.cbegin(), flow_plus.cend(), [](const RealType v)
                                           { return v >= 0.0; }));
                        // positive flow values
                        const auto flow_minus{(temp_flow.array() - temp_flow.array().abs()) / 2.0};
                        assert(flow_minus.rows() == first_coord_size + 1ll);
                        assert(std::all_of(flow_minus.cbegin(), flow_minus.cend(), [](const RealType v)
                                           { return v <= 0.0; }));

                        // const auto &flow{split_flow_field.col(j).head(first_coord_size).matrix()};
                        // exclude leftmost edge

                        A.diagonal() = 
                            A.diagonal() + 
                            flow_plus.matrix().head(first_coord_size) - 
                            flow_minus.matrix().tail(first_coord_size);
                        // exclude leftmost and rightmost edges
                        for (auto idx{1ll}; idx < A.rows(); ++idx)
                            A.coeffRef(idx, idx - 1ll) -= flow_plus(idx);
                        // exclude leftmost and rightmost edges
                        for (auto idx{0ll}; idx < A.rows() - 1ll; ++idx)
                            A.coeffRef(idx, idx + 1ll) += flow_minus(idx + 1ll);
                        // BC
                        applyBC_split_y(A, rhs, j);

                        // update current state
                        data = solve_linear_problem(A, rhs);
                    }
                }

                const SplitX splitX;
                const SplitY splitY;
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
                        const auto neibs{std::array<ptrdiff_t, 1ull>{1ll}};
                        EquationView view{A, b(0ll), 0ll, neibs};
                        bc.set_west_val(view, i);
                    }
                    {
                        const std::ptrdiff_t n{A.outerSize() - 1ll};
                        const auto neibs{std::array<ptrdiff_t, 1ull>{n - 1ll}};
                        EquationView view{A, b(n), n, neibs};
                        bc.set_east_val(view, i);
                    }
                }

                void applyBC_split_y(SpMatrix &A, RHS_t &b, ptrdiff_t j)
                {
                    {
                        const auto neibs{std::array<ptrdiff_t, 1ull>{1ll}};
                        EquationView view{A, b(0ll), 0ll, neibs};
                        bc.set_south_val(view, j);
                    }
                    {
                        const std::ptrdiff_t n{A.outerSize() - 1ll};
                        const auto neibs{std::array<ptrdiff_t, 1ull>{n - 1ll}};
                        EquationView view{A, b(n), n, neibs};
                        bc.set_north_val(view, j);
                    }
                }
            };
        } // SplittingMethod
    } // EqSolver
} // GPN