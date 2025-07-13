#pragma once
#include <cassert>
#include <iterator>

#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Properties/FaceProperties.hpp>

namespace GPN
{
    namespace Logs
    {
        namespace Rocks
        {
            struct IsPermeableLog
            {
                IsPermeableLog(
                    const auto &is_permeable_stencils,
                    const auto &grid)
                    : is_permeable{IsPermeableFactory::create(is_permeable_stencils, grid)}
                {
                    assert(is_permeable_stencils.size() == grid.mesh_size());
                }
                IsPermeable is_permeable;
            };
            struct CoreSampleLogs
            {
                CoreSampleLogs(
                    const auto &is_permeable_stencils,
                    const auto &is_perforated_stencils,
                    const auto &porosity_stencils,
                    const auto &permeability_stencils,
                    const auto &grid)
                    : is_permeable{IsPermeableFactory::create(is_permeable_stencils, grid)},
                      is_perforated{IsPerforatedFactory::create(is_perforated_stencils, is_permeable_stencils, grid)},
                      permeability{PermeabilityFactory::create(permeability_stencils, is_permeable_stencils, grid)},
                      porosity{PorosityFactory::create(porosity_stencils, is_permeable_stencils, grid)}
                {
                    assert(is_permeable_stencils.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(porosity_stencils.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(permeability_stencils.size() == grid.dual_stencils.dual_nodes.size() - 1ll);

                    assert(is_permeable.size() == grid.dual_nodes.size() - 1ll);
                    assert(permeability.size() == grid.dual_nodes.size() - 1ll);
                    assert(porosity.size() == grid.dual_nodes.size() - 1ll);

                    for (auto i{0ll}; i < is_permeable.size(); ++i)
                    {
                        assert(
                            (is_permeable(i) == 1.0) ||
                            ((is_permeable(i) == 0.0) && (is_perforated(i) == 0.0) && (porosity(i) == 0.0) && (permeability(i) == 0.0)));
                    }
                }

                IsPermeable is_permeable;
                IsPerforated is_perforated;
                Permeability permeability;
                Porosity porosity;
            };

            struct HeatLogs
            {
                HeatLogs(
                    const auto &solid_density,
                    const auto &solid_specific_heatcapacity,
                    const auto &solid_heat_conductivity,
                    const auto &porosity,
                    const auto &fluid,
                    const auto &grid)
                    : solid_density{
                          SolidDensityFactory::create(solid_density, grid)},
                      solid_specific_heatcapacity{SolidSpecificHeatCapacityFactory::create(solid_specific_heatcapacity, grid)}, 
                      medium_heat_conductivity{HeatConductivityFactory::create(porosity, solid_heat_conductivity, fluid, grid)}, 
                      solid_heat_conductivity{HeatConductivityFactory::create(solid_heat_conductivity, grid)},
                      solid_vol_heatcapacity{SolidVolumetricHeatCapacityFactory::create(solid_density, solid_specific_heatcapacity, grid)}, 
                      medium_vol_heatcapacity{MediumHeatVolumetricCapacityFactory::create(porosity, solid_density, solid_specific_heatcapacity, fluid, grid)}
                {
                    assert(solid_density.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(solid_specific_heatcapacity.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(solid_heat_conductivity.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(porosity.size() == grid.dual_stencils.dual_nodes.size() - 1ll);

                    assert(this->solid_density.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->solid_specific_heatcapacity.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->medium_heat_conductivity.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->solid_vol_heatcapacity.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->medium_vol_heatcapacity.size() == grid.dual_nodes.size() - 1ll);
                }

                const SolidDensity solid_density;
                const SolidSpecificHeatCapacity solid_specific_heatcapacity;
                const HeatConductivity solid_heat_conductivity;
                const SolidVolumetricHeatCapacity solid_vol_heatcapacity;

                const HeatConductivity medium_heat_conductivity;
                const MediumHeatVolumetricCapacity medium_vol_heatcapacity;
            };
        } // Rocks

        namespace Hydrodynamics
        {
            struct Hydrodynamics
            {
                Hydrodynamics(
                    const auto &is_permeable_stencils,
                    const auto &ext_pressure,
                    const auto &skin,
                    const auto &grid)
                    : skin{SkinFactory::create(skin, is_permeable_stencils, grid)},
                      ext_pressure{ExtPressureFactory::create(ext_pressure, is_permeable_stencils, grid)}
                {
                    assert(is_permeable_stencils.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(ext_pressure.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(skin.size() == grid.dual_stencils.dual_nodes.size() - 1ll);

                    assert(this->ext_pressure.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->skin.size() == grid.dual_nodes.size() - 1ll);

                    const auto is_permeable = IsPermeableFactory::create(is_permeable_stencils, grid);

                    for (auto i{0ll}; i < is_permeable.size(); ++i)
                    {
                        assert(
                            (is_permeable(i) == 1.0) ||
                            ((is_permeable(i) == 0.0) && (this->ext_pressure(i) == 0.0) && (this->skin(i) == 0.0)));
                    }
                }

                const ExternalPressure ext_pressure;
                const SkinFactor skin;
            };

        } // Hydrodynamics
    } // Logs

    namespace Properties
    {
        namespace Rocks
        {
            template <typename Grid2D_t>
            struct Rocks
            {
                Rocks(
                    const Logs::Rocks::CoreSampleLogs &logs,
                    const cptr<Grid2D_t> grid2D)
                    : permeability{
                          FieldFactory::create(logs.permeability, grid2D)},
                      porosity{FieldFactory::create(logs.porosity, grid2D)}
                {
                    assert(grid2D->first_coord.dual_stencils.dual_nodes.size() <= grid2D->first_coord.dual_size());
                    assert(grid2D->second_coord.dual_stencils.dual_nodes.size() <= grid2D->second_coord.dual_size());

                    assert(permeability.rows() == grid2D->first_coord.mesh_size());
                    assert(permeability.cols() == grid2D->second_coord.mesh_size());
                    assert(porosity.rows() == grid2D->first_coord.mesh_size());
                    assert(porosity.cols() == grid2D->second_coord.mesh_size());
                }

                Permeability<Grid2D_t> permeability;
                Porosity<Grid2D_t> porosity;
            };

            template <typename Grid2D_t>
            struct HeatProps
            {
                HeatProps(
                    const Logs::Rocks::HeatLogs &logs,
                    const cptr<Grid2D_t> grid2D)
                    : HeatProps{
                          logs.medium_vol_heatcapacity,
                          logs.medium_heat_conductivity,
                          grid2D}
                {
                }

                HeatProps(
                    const Logs::MediumHeatVolumetricCapacity &medium_vol_heatcapacity,
                    const Logs::HeatConductivity &heat_conductivity,
                    const cptr<Grid2D_t> grid2D)
                    : HeatProps{
                          FieldFactory::create(medium_vol_heatcapacity, grid2D),
                          FieldFactory::create(heat_conductivity, grid2D),
                          grid2D}
                {
                }

                HeatProps(
                    const MediumHeatVolumetricCapacity<Grid2D_t> &a_medium_vol_heatcapacity,
                    const MediumHeatConductivity<Grid2D_t> &a_medium_heat_conductivity,
                    const cptr<Grid2D_t> grid2D)
                    : medium_vol_heatcapacity{a_medium_vol_heatcapacity},
                      medium_heat_conductivity_axes1{a_medium_heat_conductivity},
                      medium_heat_conductivity_axes2{a_medium_heat_conductivity},
                      grid2D{grid2D}
                {
                }

                template <
                    typename Completion_t,
                    typename Well_t>
                void apply_well(
                    const Completion_t &completion,
                    const Well_t &well)
                {
#pragma region SET-HEAT-CAPACITY
                    // first column -- inside the tube, contains only water
                    medium_vol_heatcapacity.col(0ll) /*.head(tube_end)*/ =
                        completion.flow().volumetric_heat_capacity;
                    // second column -- from tube inner radius to column outer radius
                    medium_vol_heatcapacity.col(1ll) /*.head(tube_end)*/ =
                        completion.casing_volumetric_heat_capacity();
                    // third column -- cement cross-section
                    medium_vol_heatcapacity.col(2ll) /*.head(tube_end)*/ =
                        completion.cement_volumetric_heat_capacity();
#pragma endregion
#pragma region SET-HEAT-CONDUCTIVITY
                    // heat conductivity of flowing water in r-direction is infinite
                    this->medium_heat_conductivity_axes2.col(0ll) = 
                        std::numeric_limits<RealType>::infinity();
                    // put values for cementOuter at medium_vol_heatcapacity.col(1ll).
                    // CementOuter is a part of col(1ll)
                    const auto &grid_r = grid2D->second_coord;
                    const auto &sandface = completion.back();
                    const ptrdiff_t id{1ll};
                    medium_heat_conductivity_axes2.col(1ll) =
                        sandface.heat_conductivity /
                        std::log(sandface.outer_radius / sandface.inner_radius) *
                        std::log(grid_r.dual_nodes(id + 1ll) / grid_r.mesh_nodes(id));
                    // r_{1/2} is fixed at HeatFaceProps container

                    // interpolate verticle heat conductivity:
                    // (1) modify water heat conductivity in col(0ll)
                    const auto &flow = completion.front();
                    medium_heat_conductivity_axes1.col(0ll) =
                        flow.heat_conductivity * flow.area() / grid2D->face_area_axes1(0ll);
                    for (const auto v : grid2D->face_area_axes1)
                        assert(flow.area() <= v + 1e-12);
                    // (2) set sandwich heat conductivity in col(1ll)
                    medium_heat_conductivity_axes1.col(1ll) =
                        completion.integral_vertical_casing_heat_conductivity();
                    // (3) set cement heat conductivity in col(2ll)
                    medium_heat_conductivity_axes1.col(2ll) =
                        completion.integral_vertical_cement_heat_conductivity();
#pragma endregion
                }

                MediumHeatVolumetricCapacity<Grid2D_t> medium_vol_heatcapacity;
                // orthotropic medium
                MediumHeatConductivity<Grid2D_t>
                    medium_heat_conductivity_axes1,
                    medium_heat_conductivity_axes2;
                const cptr<Grid2D_t> grid2D;
            };
        } // Rocks
    } // Properties

    namespace FaceProperties
    {
        namespace Rocks
        {
            template <typename Grid2D_t>
            struct HeatFaceProps
            {
                HeatFaceProps(
                    const Properties::Rocks::HeatProps<Grid2D_t> &props,
                    const cptr<Grid2D_t> grid2D)
                    : medium_heat_conductivity{
                          FaceInterpolatedFieldFactory::create(
                              props.medium_heat_conductivity_axes1,
                              props.medium_heat_conductivity_axes2,
                              grid2D)},
                      grid2D{grid2D}
                {
                }

                template <
                    typename Completion_t,
                    typename Well_t>
                void apply_well(
                    const Completion_t &completion,
                    const Well_t &well)
                {
#pragma region SET-HEAT-CONDUCTIVITY
                    medium_heat_conductivity.face_vals_axes2.col(0ll) /*.head(tube_end)*/ =
                        completion.integral_inner_radial_heat_conductivity();
#pragma endregion
                }

                MediumHeatConductivity<Grid2D_t> medium_heat_conductivity;
                const cptr<Grid2D_t> grid2D;
            };
        } // Rocks

    } // FaceProperties
} // GPN