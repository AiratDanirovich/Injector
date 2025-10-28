#pragma once

#include <nlohmann/json.hpp>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

using json = nlohmann::json;

GPN::PhasePropertiesJT make_water(const json &data)
{
  using namespace GPN;
  using namespace GPN::Phases;

  /*fluid*/
  RealType
      viscosity{data["fluid"]["viscosity"]},
      density{data["fluid"]["density"]},
      capacity{data["fluid"]["specific_heat_capacity"]},
      heat_conductivity{data["fluid"]["heat_conductivity"]},
      joule_thomson{data["fluid"]["joule_thomson"]};

  return PhasePropertiesJT{
      FluidFactory::create_water_JT(
          Viscosity{viscosity},
          GPN::Density{density},
          GPN::SpecificHeatCapacity{capacity},
          GPN::HeatConductivity{heat_conductivity},
          JouleThomson{joule_thomson})};
}