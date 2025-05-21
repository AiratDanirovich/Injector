#include <iostream>
#include <memory>

#include <Injector/Grids/Defines.h>
#include <Injector/Solver/SolverFactory.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

TEST_CASE("Solver")
{
    RealType val{1.0};
    SolverFactory solver_factory{val, Box{Segment{0, 1}, Segment{0, 1}}, 5, 11};
    const auto grid2D{solver_factory.grid()};
    const auto conductivity_field{solver_factory.conductivity_field()};
    const auto capacity_field{solver_factory.capacity_field()};
    const auto initial_state{solver_factory.initial_state()};
    const auto bc{solver_factory.boundary_conditions()};

    const RealType well_rate{1.0};

    const auto flow_field{solver_factory.flow_field(well_rate)};

    Solver solver{
        conductivity_field,
        flow_field, grid2D,
        capacity_field,
        initial_state,
        bc, 0.0};

    solver.advance(0.005);
}