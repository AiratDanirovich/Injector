// #include <memory>
#include <iostream>
#include <fstream>
#include <vector>
// #include <string>
// #include <numbers>
// #include <cmath>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/Grids2D.hpp>
// #include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Model/WellFactory.hpp>
// #include <Injector/Properties/Factory.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace Catch;
using namespace Catch::Matchers;
using json = nlohmann::json;

using namespace std;
using namespace GPN;
using namespace GPN::Logs;
using namespace GPN::Grids;
using namespace GPN::Phases;
using namespace GPN::Completion;

using VR = std::vector<RealType>;

LogValuesContainer transfer_to_eigen(const VR &data, const RealType factor = 1.0)
{
  LogValuesContainer out(data.size());
  for (auto i{0ull}; i < data.size(); ++i)
    out(i) = factor * data[i];
  return out;
}

// VR generate_stencils(RealType t0, RealType t1, RealType t_step_major)
// {
//   auto segm_count{static_cast<size_t>(std::ceil(t1 - t0) / t_step_major)};
//   double step = (t1 - t0) / (segm_count);
//   VR out(segm_count + 1ll);

//   for (auto i{0ull}; i < out.size(); ++i)
//     out[i] = t0 + i * step;
//   return out;
// }

// VR generate_steps(const VR &dual_nodes)
// {
//   VR out(dual_nodes.size() - 1ll);

//   for (auto i{0ull}; i < out.size(); ++i)
//     out[i] = dual_nodes[i + 1] - dual_nodes[i];
//   return out;
// }

const VR make_r_stencils(const json &data, const auto &well_holes)
{
  const std::string r_grid_type = data["grid"]["r_grid_type"];
  const RealType
      rMin{data["grid"]["r_start"]},
      rMax{data["grid"]["r_end"]};

  if (r_grid_type == "uniform")
  {
    const auto &data2 = data["grid"]["r_uniform_grid"];
    return well_holes.generate_uniform_radial_grid(
        rMin, rMax, data2["rNodes"]);
  }
  else if (r_grid_type == "log")
  {
    const auto &data2 = data["grid"]["r_log_grid"];
    return well_holes.generate_log_radial_grid(
        rMin, rMax, data2["q"], data2["r_max_step"]);
  }
  else
    throw std::runtime_error("Incorrect radial grid descriptor.");
}


const auto get_completion(const json &data)
{
  using namespace GPN::Completion;

  std::vector<Ring> out;
  out.reserve(6);

  { // flowing fluid
    const auto &data2 = data["fluid"];
    out.push_back(
        Ring{
            Flow{
                Density{data2["density"]},
                SpecificHeatCapacity{data2["specific_heat_capacity"]},
                GPN::HeatConductivity{data2["heat_conductivity"]}},
            Thickness{data["completion"]["tube"]["inner_radius"]},
            InnerRadius{0.0},
            Depth{std::numeric_limits<RealType>::max()}});
  }

  { // tube
    const auto &data2 = data["completion"]["tube"];

    out.push_back(
        Ring{
            Tube{
                Density{data2["density"]},
                SpecificHeatCapacity{data2["specific_heat_capacity"]},
                GPN::HeatConductivity{data2["heat_conductivity"]}},
            Thickness{data2["thickness"]},
            InnerRadius{out.back().outer_radius},
            Depth{data2["depth"]}});
  }
  
  { // annulus
    const auto &data2 = data["completion"]["annulus"];

    out.push_back(
        Ring{
            Annulus{
                Density{data2["density"]},
                SpecificHeatCapacity{data2["specific_heat_capacity"]},
                GPN::HeatConductivity{data2["heat_conductivity"]}},
            Thickness{data2["thickness"]},
            InnerRadius{out.back().outer_radius},
            Depth{data2["depth"]}});
  }
    
  { // column
    const auto &data2 = data["completion"]["column"];

    out.push_back(
        Ring{
            Annulus{
                Density{data2["density"]},
                SpecificHeatCapacity{data2["specific_heat_capacity"]},
                GPN::HeatConductivity{data2["heat_conductivity"]}},
            Thickness{data2["thickness"]},
            InnerRadius{out.back().outer_radius},
            Depth{data2["depth"]}});
  }
      
  { // cement_1
    const auto &data2 = data["completion"]["cement"]["inner"];

    out.push_back(
        Ring{
            Cement{
                Density{data2["density"]},
                SpecificHeatCapacity{data2["specific_heat_capacity"]},
                GPN::HeatConductivity{data2["heat_conductivity"]}},
            Thickness{data2["thickness"]},
            InnerRadius{out.back().outer_radius},
            Depth{data2["depth"]}});
  }
      
  { // cement_2
    const auto &data2 = data["completion"]["cement"]["outer"];

    out.push_back(
        Ring{
            Cement{
                Density{data2["density"]},
                SpecificHeatCapacity{data2["specific_heat_capacity"]},
                GPN::HeatConductivity{data2["heat_conductivity"]}},
            Thickness{data2["thickness"]},
            InnerRadius{out.back().outer_radius},
            Depth{data2["depth"]}});
  }

  return out;
}

TEST_CASE("Solver", "SelfSimilarCyl")
{
  ifstream f("heatflow_test_data.json");
  REQUIRE(f.is_open());
  json data = json::parse(f);

  /*START*/
  // input parameters
  /*fluid*/
  RealType
      viscosity{data["fluid"]["viscosity"]},
      density{data["fluid"]["density"]},
      capacity{data["fluid"]["specific_heat_capacity"]},
      heat_conductivity{data["fluid"]["heat_conductivity"]};
  /*collector*/
  const VR thickness = data["collector"]["thickness"];
  // hydrodynamic logs
  const auto is_permeable_stencils{transfer_to_eigen(data["collector"]["is_permeable"])};
  const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"])};
  const auto porosity_stencils{transfer_to_eigen(data["collector"]["porosity"])};
  const auto permeability_stencils{transfer_to_eigen(data["collector"]["permeability"], 1e-12)};
  const auto weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"])};
  // heat logs
  const auto solid_heatconductivity_stencils{transfer_to_eigen(data["collector"]["heatConductivity"])};
  const auto solid_density_stencils{transfer_to_eigen(data["collector"]["solidDensity"])};
  const auto solid_specific_heatcapacity_stencils{transfer_to_eigen(data["collector"]["solidSpecificHeatCapacity"])};
  /*grid*/
  const RealType
      z_minor_step{data["grid"]["z_minor_step"]}; // m
  /*completion*/
  const Casing completion{get_completion(data)};
  /*END*/

  // make grid2D
  // r_stencils
  const WellHoles well_holes{WellHolesFactory::create(completion)};
  const VR r_stencils{make_r_stencils(data, well_holes)};
  const RealType &rMax = r_stencils.back();
  const RealType &rMin = r_stencils.front();

  // z-refiner
  RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
  // the grid itself
  const auto grid2D{
      Grids::CylinderGridFactory::create(refiner,
                                         Grids::Factory::generate_dual_grid_stencils_from_steps(
                                             0.0, thickness),
                                         r_stencils)};
  const auto &grid{grid2D->first_coord};

  cout << "radial dual grid stencils:\n"
       << transfer_to_eigen(r_stencils).transpose() << endl;

  const Logs::Rocks::CoreSampleLogs core_data{
      is_permeable_stencils,
      is_perforated_stencils,
      porosity_stencils,
      permeability_stencils,
      grid};

  // make fluid
  const PhaseProperties water{
      FluidFactory::create_water(
          Viscosity{viscosity},
          Density{density},
          SpecificHeatCapacity{capacity},
          GPN::HeatConductivity{heat_conductivity})};
  const auto weights{
      RateWeightsFactory::create(
          weights_stencils,
          core_data.is_permeable,
          grid)};
  const Well_Explicit well{
      core_data.is_permeable, core_data.is_perforated, weights};

  const Logs::Rocks::HeatLogs heat_logs{
      solid_density_stencils,
      solid_specific_heatcapacity_stencils,
      solid_heatconductivity_stencils,
      porosity_stencils,
      water,
      grid2D->first_coord};

  Properties::Rocks::HeatProps heat_props{
      heat_logs, grid2D};

  // properties of material that fills the well up to the sandface
  heat_props.apply_well(completion, well);

  FaceProperties::Rocks::HeatFaceProps heat_face_props{
      heat_props, grid2D};
  heat_face_props.apply_well(completion, well);
}