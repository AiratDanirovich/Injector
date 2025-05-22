#include <Injector/Grids/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Solver/SplittingMethod/BaseSplit.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::EqSolver::SplittingMethod;

using VR = std::vector<RealType>;

// input data
const VR z_stencils = Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5);
const VR r_stencils = Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11);
const VR conductivity_stencils = Logs::StencilsFactory::generate_conductivity_StepProperty(z_stencils);

// Tests Cylinder grid, (r; z)
TEST_CASE("BaseSplitTest")
{
    const Grids::CylinderGridFactory grid_factory{z_stencils, r_stencils};
    const Logs::HeatLogsFactory heat_factory{conductivity_stencils, grid_factory.grid()->first_coord};

    Properties::HeatConductivity conductivity_field{
        heat_factory.conductivity,
        grid_factory.grid()};

    BaseSplit base_split{
        conductivity_field.face_vals_axes2,
        grid_factory.grid()->first_coord.mesh_size(),
        grid_factory.grid()->second_coord.mesh_size()};
}
