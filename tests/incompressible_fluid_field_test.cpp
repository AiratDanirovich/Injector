
#include <vector>
#include <fstream>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/Grids2D.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well/WellHoles.hpp>
#include <Injector/Model/Well/WellFactory.hpp>
#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Hydrodynamic/Incompressible/IncompressibleFluid.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/set_is_permeable_stencils.hpp"
#include "includes/make_r_stencils.hpp"
#include "includes/get_completion.hpp"

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;
using json = nlohmann::json;

using namespace std;
using namespace GPN;
using namespace GPN::Logs;
using namespace GPN::CrossFlow;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::Completion;
using namespace GPN::Hydrodynamic;

using VR = std::vector<GPN::RealType>;

TEST_CASE("Solver", "SelfSimilarCyl")
{

  ifstream f("heatflow_test_data.json");
  REQUIRE(f.is_open());
  json data = json::parse(f);

  /*fluid*/
  RealType
      viscosity{data["fluid"]["viscosity"]},
      density{data["fluid"]["density"]},
      capacity{data["fluid"]["specific_heat_capacity"]},
      heat_conductivity{data["fluid"]["heat_conductivity"]},
      joule_thomson{data["fluid"]["joule_thomson"]};

  /*collector*/
  const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
  const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};

  const RealType zTop{0.0};
  const auto
      z_minor_step{data["grid"]["z_minor_step"].get<RealType>()}; // m
  const auto thickness{data["collector"]["thickness"].get<VR>()};
  const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"].get<VR>(), 1e-12)};
  const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
  const auto is_permeable_stencils{set_is_permeable_stencils(is_perforated_stencils, to_layers)};

    /*completion*/
  // z-refiner
  const auto z_stencils{Grids::Factory::generate_dual_grid_stencils_from_steps(
      0.0, thickness)};
  RefinerVerticle z_refiner{z_minor_step, is_permeable_stencils};
  const auto temp_grid_z{
      Factory::create_axes<CoordinateTypes::Z>(
          z_refiner, z_stencils)};
  const Casing<VarRing> completion{get_completion(data, temp_grid_z)};

  // make grid2D
  // r_stencils
  const VR r_stencils{
      WellHoles{WellHolesFactory::create(completion)}.get_stencils(
          data["grid"]["r_start"],
          data["grid"]["r_end"])};
  // r-refiner
  const AbstractRefinerRadial *r_refiner{
      make_r_refiner(data)};

  const auto r_nodes{r_refiner->refine(r_stencils)};
  // the grid itself
  const auto grid2D{
      Grids::CylinderGridFactory::create(
          z_refiner, z_stencils,
          r_nodes)};
  const auto &grid_z{grid2D->first_coord};
  const auto &grid_r{grid2D->second_coord};

  const auto rMin{grid_r.dual_front()};
  const auto rMax{grid_r.dual_back()};

  const auto is_permeable{
      Logs::IsPermeableFactory::create(
          is_permeable_stencils,
          grid_z)};

  const auto is_perforated{
      Logs::IsPerforatedFactory::create(
          is_perforated_stencils,
          is_permeable_stencils,
          grid_z)};

  const auto permeability{
      Logs::PermeabilityFactory::create(
          permeability_stencils,
          is_permeable_stencils,
          grid_z)};

  // make fluid
  const PhasePropertiesJT water{
      FluidFactory::create_water_JT(
          Viscosity{viscosity},
          GPN::Density{density},
          GPN::SpecificHeatCapacity{capacity},
          GPN::HeatConductivity{heat_conductivity},
          JouleThomson{joule_thomson})};

  const auto pressure_field{IncompressibleFluidField{
      water,
      permeability,
      grid2D
  }};
}