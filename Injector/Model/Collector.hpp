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
                      solid_specific_heatcapacity{SolidSpecificHeatCapacityFactory::create(solid_specific_heatcapacity, grid)}, medium_heat_conductivity{HeatConductivityFactory::create(porosity, solid_heat_conductivity, fluid, grid)}, solid_vol_heatcapacity{SolidVolumetricHeatCapacityFactory::create(solid_density, solid_specific_heatcapacity, grid)}, medium_vol_heatcapacity{MediumHeatVolumetricCapacityFactory::create(porosity, solid_density, solid_specific_heatcapacity, fluid, grid)}
                {
                    assert(solid_density.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(solid_specific_heatcapacity.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(heat_conductivity.size() == grid.dual_stencils.dual_nodes.size() - 1ll);
                    assert(porosity.size() == grid.dual_stencils.dual_nodes.size() - 1ll);

                    assert(this->solid_density.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->solid_specific_heatcapacity.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->heat_conductivity.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->solid_vol_heatcapacity.size() == grid.dual_nodes.size() - 1ll);
                    assert(this->medium_vol_heatcapacity.size() == grid.dual_nodes.size() - 1ll);
                }

                const SolidDensity solid_density;
                const SolidSpecificHeatCapacity solid_specific_heatcapacity;
                const HeatConductivity medium_heat_conductivity;
                const SolidVolumetricHeatCapacity solid_vol_heatcapacity;
                const MediumHeatVolumetricCapacity medium_vol_heatcapacity;

                // private:
                //     template <typename T>
                //     static T multiply(
                //         const T &lhs,
                //         const T &rhs)
                //     {
                //         T out(lhs.size(), 0.0);
                //         for (auto i{0ll}; i < lhs.size(); ++i)
                //             out[i] = lhs[i] * rhs[i];

                //         return out;
                //     }
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
                    const Logs::MediumHeatVolumetricCapacity &medium_vol_heatcapacity,
                    const Logs::HeatConductivity &heat_conductivity,
                    const cptr<Grid2D_t> grid2D)
                    : medium_vol_heatcapacity{
                          FieldFactory::create(medium_vol_heatcapacity, grid2D)},
                      medium_heat_conductivity{FieldFactory::create(heat_conductivity, grid2D)}
                {
                }

                HeatProps(
                    const Logs::Rocks::HeatLogs &logs,
                    const cptr<Grid2D_t> grid2D)
                    : medium_vol_heatcapacity{
                          FieldFactory::create(logs.medium_vol_heatcapacity, grid2D)},
                      medium_heat_conductivity{FieldFactory::create(logs.medium_heat_conductivity, grid2D)}
                {
                }

                template <
                    typename Completion_t,
                    typename Well_t,
                    typename Fluid_t>
                void apply_well(
                    const Completion_t &completion,
                    const Well_t &well,
                    const Fluid_t &fluid)
                {
#pragma region SET-HEAT-CAPACITY
                    // first column -- inside the tube, contains only water
                    medium_vol_heatcapacity.col(0ll) /*.head(tube_end)*/ =
                        fluid.volumetric_heat_capacity;
                    // second column -- from tune inner radius to sandface radius
                    medium_vol_heatcapacity.col(1ll) /*.head(tube_end)*/ =
                        completion.volumetric_heat_capacity();

                    // medium_vol_heatcapacity.col(0ll).tail(medium_vol_heatcapacity.rows() - tube_end) =
                    //     fluid.volumetric_heat_capacity*(r_column*r_column)/(r_tube*r_tube);

                    // medium_vol_heatcapacity.col(1ll).head(tube_end) =
                    //     upper_annulus_capacity;
                    // medium_vol_heatcapacity.col(1ll).tail(medium_vol_heatcapacity.rows() - tube_end) =
                    //     lower_annulus_capacity;
#pragma endregion
#pragma region SET-HEAT-CONDUCTIVITY
                    // put values for cementOuter at medium_vol_heatcapacity.col(1ll).
                    // CementOuter is a part of col(1ll)
                    const auto &grid = grid2D->first_coord;
                    const auto &sandface = completion.back();
                    const ptrdiff_t id{1ll};
                    medium_heat_conductivity.col(1ll) =
                        sandface.heat_conductivity /
                        std::log(sandface.outer_radius / sandface.inner_radius) *
                        std::log(grid.dual_nodes(id + 1ll) / grid.mesh_nodes(id));
                    // r_{1/2} is fixed at HeatFaceProps container
#pragma endregion
                }

                MediumHeatConductivity<Grid2D_t> medium_heat_conductivity;
                MediumHeatVolumetricCapacity<Grid2D_t> medium_vol_heatcapacity;
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
                              props.medium_heat_conductivity,
                              grid2D)},
                      grid2D{grid2D}, props{props}
                {
                }

                template <
                    typename Completion_t,
                    typename Well_t,
                    typename Fluid_t>
                void apply_well(
                    const Completion_t &completion,
                    const Well_t &well,
                    const Fluid_t &fluid)
                {
                    // last row with the tube
                    const auto &mesh = grid2D->first_coord.mesh_nodes;
                    // const auto it = std::upper_bound(mesh.cbegin(), mesh.cend(), completion.tube_depth);
                    // const ptrdiff_t tube_end{std::distance(mesh.cbegin(), it) - 1ll};

#pragma region SET-HEAT-CONDUCTIVITY
                    medium_heat_conductivity.face_vals_axes2.col(0ll) /*.head(tube_end)*/ =
                        completion.integral_inner_radial_heat_conductivity();
#pragma endregion
                }

                MediumHeatConductivity<Grid2D_t> medium_heat_conductivity;
                const Properties::Rocks::HeatProps<Grid2D_t> &props;
                const cptr<Grid2D_t> grid2D;
            };
        } // Rocks

    } // FaceProperties
} // GPN