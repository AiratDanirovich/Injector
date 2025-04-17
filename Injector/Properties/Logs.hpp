#pragma once

// #include <iostream>

#include <Eigen/Core>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace Logs
    {
        using StepPropertyContainer = Eigen::ArrayX<RealType>;
        using InternalFaceValues = Eigen::ArrayX<RealType>;

        struct InterpolatedDataContainer : public StepPropertyContainer
        {
            using StepPropertyContainer::StepPropertyContainer;
        };

        /// @brief Container for values of step properties.
        struct StepProperty
        {
            StepProperty(
                const std::vector<RealType> &adata)
                : data(adata.size())
            {
#pragma region ASSERTIONS
                assert(data.size() > 0ull);
                for (auto idx{0ull}; idx < adata.size(); ++idx)
                    // all properties are non-negative
                    assert(adata[idx] >= 0.0);
#pragma endregion
                for (auto idx{0ull}; idx < adata.size(); ++idx)
                    data(idx) = adata[idx];
            }

            StepPropertyContainer data;
        };

        struct StepPropertyGrid
        {
            using Grid_t = Grids::GridDual; // AxesGrid<CoordinateTypes::Z>;

            StepPropertyGrid(
                const StepProperty &step_vals,                      // stencil values per every layer
                const InterpolatedDataContainer &interpolated_vals, // interpolated values corresponding to the grid
                const Grid_t &grid)
                : property_vals{step_vals}, 
                log_vals{interpolated_vals}, 
                grid{grid}
            {
            }

            StepPropertyGrid(
                const StepProperty &property_vals,
                const Grid_t &grid)
                : StepPropertyGrid{
                      // stencil values per every layer
                      property_vals,
                      // interpolated stencil-values
                      // on refined grid nodes
                      interpolate(property_vals, grid),
                      grid}
            {
            }

            auto operator()(auto id) const
            {
                assert(id < log_vals.size());
                return log_vals(id);
            }

            const InterpolatedDataContainer &log_vals;
            const Grid_t &grid;

        protected:
            // values between stencil nodes
            StepProperty property_vals;

        private:
            static InterpolatedDataContainer interpolate(
                const StepProperty &property_vals,
                const Grid_t &grid)
            {
                InterpolatedDataContainer out(grid.dual_size());
                // interpolate property_vals on the grid
                for (
                    // loop through all control_volumes
                    auto volume_id{0ll}, mesh_node_id{0ll};
                    volume_id < grid.dual_steps.size();
                    ++volume_id)
                {
                    //    std::cout << "volume_id =     " << volume_id << std::endl;
                    //    std::cout << "mesh_node_ids = " << std::endl;
                    // set constant value within a fixed control volume
                    for (;
                         (mesh_node_id < grid.dual_steps.size()) &&
                         (grid.mesh_nodes(mesh_node_id) < grid.dual_stencils(volume_id + 1ull));
                         ++mesh_node_id)
                    {
                        //        std::cout << mesh_node_id << ' ';
                        out(mesh_node_id) = property_vals.data(volume_id);
                    }
                    //    std::cout << std::endl;
                }

                return out;
            }
        };

        template <typename CoordinateType_t>
        struct FaceInterpolator
        {
            using Axes = CoordinateType_t;

            const InternalFaceValues face_values;

            FaceInterpolator(
                const CoordinateType_t &coord,
                const StepPropertyGrid &log)
                : face_values{interpolate<CoordinateType_t>(coord, log)}
            {
            }

        protected:
            template <typename CoordinateType_t>
            static auto interpolate(
                const CoordinateType_t &coord,
                const StepPropertyGrid &log)
            {
                const auto &grid{log.grid};
                // if dual_size() == 0 --- no internal faces
                assert(grid.dual_size() > 2ll);

                InternalFaceValues out(coord.dual_size() - 2ll);

                for (auto id{0ll}; id < out.size(); ++id)
                    out(id) = CoordinateType_t::face_interpolator(
                        grid.mesh_nodes(id), grid.mesh_nodes(id + 1ll), grid.dual_nodes(id + 1ll), log(id), log(id + 1ll));

                return out;
            }
        };

        template <typename CoordinateType_t /* = CoordinateTypes::Z*/>
        struct FaceInterpolatedProperty
            : public StepPropertyGrid,
              public FaceInterpolator<CoordinateType_t>
        {
            FaceInterpolatedProperty(
                const StepPropertyGrid &vals,
                const FaceInterpolator<CoordinateType_t> &interp) noexcept
                : StepPropertyGrid{vals},
                  FaceInterpolator<CoordinateType_t>{interp}
            {
            }
        };

        /// @brief Property that must only contain {0; 1} values
        struct IndicatorProperty : public StepPropertyGrid
        {
            IndicatorProperty(const StepPropertyGrid &vals)
                : StepPropertyGrid{vals}
            {
                for (auto idx{0ll}; idx < vals.log_vals.size(); ++idx)
                    assert(log_vals(idx) == 0.0 || log_vals(idx) == 1.0);
            }
        };

        struct IsPermeable : public IndicatorProperty
        {
            using IndicatorProperty::IndicatorProperty;
        };

        struct Permeability : public StepPropertyGrid
        {
            using StepPropertyGrid::StepPropertyGrid;
        };
        struct Porosity : public StepPropertyGrid
        {
            using StepPropertyGrid::StepPropertyGrid;
        };
        struct SkinFactor : public StepPropertyGrid
        {
            using StepPropertyGrid::StepPropertyGrid;
        };
        struct HeatConductivity : public FaceInterpolatedProperty<CoordinateTypes::Z>
        {
            using FaceInterpolatedProperty<
                    CoordinateTypes::Z>
                ::FaceInterpolatedProperty;
        };
        struct MatrixHeatCapacity : public StepPropertyGrid
        {
            using StepPropertyGrid::StepPropertyGrid;
        };
        struct MatrixDensity : public StepPropertyGrid
        {
            using StepPropertyGrid::StepPropertyGrid;
        };
    } // Logs
} // GPN