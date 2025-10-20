#include <iostream>
#include <memory>

#include <Injector/Grids/Defines.h>
#include <Injector/Properties/FaceProperties.hpp>

#include <Injector/History/ZeroRatesFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/SplittingMethod/Solver.hpp>

#include "includes/BCFunctor.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace GPN;
using namespace std;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::SplittingMethod;

// const RealType well_rate{1.0};
const RealType val{1.0};
const auto z_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 5)};
const auto r_stencils{
    Grids::Factory::generate_dual_grid_stencils_uniform(0, 1, 11)};

// hydrodynamic logs
const auto is_permeable_stencils{
    Logs::RawDataFactory::generate_is_permeable(z_stencils)};
auto is_perforated_stencils{
    Logs::RawDataFactory::generate_is_permeable(z_stencils)};
const auto porosity_stencils{
    Logs::RawDataFactory::generate_porosity(z_stencils, is_permeable_stencils)};
const auto permeability_stencils{
    Logs::RawDataFactory::generate_permeability(z_stencils, is_permeable_stencils)};
const auto ext_pressure_stencils{
    Logs::RawDataFactory::generate_ext_pressure(z_stencils, is_permeable_stencils)};
// heat logs
const auto solid_density_stencils{
    Logs::RawDataFactory::generate_solid_density(z_stencils)};
const auto solid_specific_heatcapacity_stencils{
    Logs::RawDataFactory::generate_solid_specific_heatcapacity(z_stencils)};
const auto solid_heatconductivity_stencils{
    Logs::RawDataFactory::generate_conductivity(z_stencils)};

TEST_CASE("Solver")
{
    auto it = std::ranges::find_if(
        is_perforated_stencils,
        [](RealType v)
        { return v == 1.0; });
    (*it) = 0.0;

    const auto grid2D{Grids::CylinderGridFactory::create(z_stencils, r_stencils)};
    const auto &grid_z{grid2D->first_coord()};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        solid_heatconductivity_stencils,
        porosity_stencils,
        Phases::FluidFactory::create_water(1.0, 1.0),
        grid_z};

    const Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};

    const FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};

    const Problem::InitialCondition initial_state{
        State::State2D::FillWithConst(
            *grid2D, val)};

    const BoundaryConditions::GeneralBC bc{
        grid2D,
        make_shared<BCFunctor>(val),
        BoundaryConditions::GeneralBC::BoundaryCondition::first};

    const Logs::Rocks::CoreSampleLogs core_data{
        is_permeable_stencils,
        is_perforated_stencils,
        porosity_stencils,
        permeability_stencils,
        grid_z};

  // external pressure log
  const auto external_pressure{
      Logs::ExtPressureFactory::create(
          ext_pressure_stencils,
          is_permeable_stencils,
          grid_z)};
  // rates field factory
  FaceProperties::ZeroRatesFactory rates_factory{
      grid2D, core_data.is_permeable, external_pressure};

    Solver solver{
        heat_face_props.medium_heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        rates_factory, initial_state,
        bc, 0.0};

    solver.advance(0.005);
}