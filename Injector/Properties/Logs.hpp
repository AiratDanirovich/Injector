#pragma once

// #include <iostream>

#include <Eigen/Core>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>
#include <Injector/Grids/Grids1D.hpp>

#include <Injector/Model/Phases/PhaseProperties.hpp>

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
            StepPropertyContainer data;

            StepProperty(
                const StepPropertyContainer &adata)
                : StepProperty{
                      std::vector<RealType>(
                          adata.begin(),
                          adata.end())}
            {
            }

            StepProperty(
                const std::vector<RealType> &adata)
                : data(adata.size())
            //: data(adata.begin(), adata.end())
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

            auto size() const
            {
                return data.size();
            }

            auto operator()(auto id) const
            {
                return data(id);
            }

            auto operator*(RealType c) const{
                const auto &lhs{*this};
                return StepProperty{lhs.data*c};

            }
            auto operator*(const StepProperty &rhs) const
            {
                const auto &lhs{*this};
                return StepProperty{lhs.data * rhs.data};
            }

            auto operator+(const StepProperty &rhs) const
            {
                const auto &lhs{*this};
                return StepProperty{lhs.data + rhs.data};
            }
        };

        
        auto operator-(RealType v, const StepProperty &rhs)
        {
            return StepProperty{v - rhs.data};
        }



        struct AssertNonNegative
        {
            AssertNonNegative(const StepProperty &vals)
            {
                for (std::ptrdiff_t id{0ll}; id < vals.size(); ++id)
                    assert(vals(id) >= 0.0);
            }
        };

        struct StepPropertyGrid
        {
            using Grid_t = Grids::GridDual; // AxesGrid<CoordinateTypes::Z>;

            StepPropertyGrid(const StepPropertyGrid &) noexcept = default;
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

            auto stencil_size() const
            {
                return property_vals.size();
            }

            auto size() const
            {
                return log_vals.size();
            }

            operator const StepProperty &() const { return property_vals; }

            const InterpolatedDataContainer log_vals;
            const Grid_t &grid;
            // values between stencil nodes
            const StepProperty property_vals;

        private:
            static InterpolatedDataContainer interpolate(
                const StepProperty &property_vals,
                const Grid_t &grid)
            {
                InterpolatedDataContainer out(grid.mesh_size());
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

        /// @brief Property that must only contain {0; 1} values
        struct IndicatorProperty
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            IndicatorProperty(
                const StepPropertyGrid &vals)
                : StepPropertyGrid{vals},
                  AssertNonNegative{vals}
            {
                for (auto idx{0ll}; idx < vals.log_vals.size(); ++idx)
                    assert(log_vals(idx) == 0.0 || log_vals(idx) == 1.0);
            }
        };

        struct IsPermeable : public IndicatorProperty
        {
            using IndicatorProperty::IndicatorProperty;
        };

        struct Permeability
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            Permeability(
                const StepPropertyGrid &permeability,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{permeability},
                  AssertNonNegative{permeability}
            {
                assert(permeability.size() == is_permeable.size());
                for (std::ptrdiff_t id{0ll}; id < permeability.size(); ++id)
                    assert(
                        ((is_permeable(id) == 1.0) && (permeability(id) > 0.0)) ||
                        ((is_permeable(id) == 0.0) && (permeability(id) == 0.0)));
            }
        };

        struct Porosity
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            Porosity(
                const StepPropertyGrid &porosity,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{porosity},
                  AssertNonNegative{porosity}
            {
                assert(porosity.size() == is_permeable.size());
                for (std::ptrdiff_t id{0ll}; id < porosity.size(); ++id)
                {
                    assert(
                        ((is_permeable(id) == 1.0) && (porosity(id) > 0.0)) ||
                        ((is_permeable(id) == 0.0) && (porosity(id) == 0.0)));
                    assert(porosity(id) <= 1.0);
                }
            }
        };
        struct SkinFactor : public StepPropertyGrid
        {
            SkinFactor(
                const StepPropertyGrid &skin,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{skin}
            {
                assert(skin.size() == is_permeable.size());
                for (std::ptrdiff_t id{0ll}; id < skin.size(); ++id)
                    assert(
                        ((is_permeable(id) == 1.0)) ||
                        ((is_permeable(id) == 0.0) && (skin(id) == 0.0)));
            }
        };

        struct SolidDensity
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            SolidDensity(
                const StepPropertyGrid &density,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{density},
                  AssertNonNegative{density}
            {
                assert(density.size() == is_permeable.size());
            }
        };
        struct SolidSpecificHeatCapacity
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            SolidSpecificHeatCapacity(
                const StepPropertyGrid &capacity,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{capacity},
                  AssertNonNegative{capacity}
            {
                assert(capacity.size() == is_permeable.size());
            }
        };

        struct SolidVolumetricHeatCapacity
            : public StepPropertyGrid
        {
            static_assert(
                std::is_same<
                    SolidDensity::Grid_t,
                    SolidSpecificHeatCapacity::Grid_t>::value);

            SolidVolumetricHeatCapacity(
                const SolidDensity &density,
                const SolidSpecificHeatCapacity &heat_capacity)
                : StepPropertyGrid{
                      static_cast<StepProperty>(density) * static_cast<StepProperty>(heat_capacity),
                      density.grid}
            {
            }
        };

        template <typename CoordinateType_t /* = CoordinateTypes::Z*/>
        struct FaceInterpolator
        {
            using Axes = CoordinateType_t;

            const InternalFaceValues face_values;

            FaceInterpolator(
                const StepPropertyGrid &log)
                : face_values{
                      interpolate(log)}
            {
            }

        protected:
            static auto interpolate(
                const StepPropertyGrid &log)
            {
                const auto &grid{log.grid};
                // if dual_size() == 0 --- no internal faces
                assert(grid.dual_size() > 2ll);

                InternalFaceValues out(log.grid.dual_size() - 2ll);

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
                const StepPropertyGrid &vals) noexcept
                : StepPropertyGrid{vals},
                  FaceInterpolator<CoordinateType_t>{vals}
            {
            }
        };

        using ZInterpolator =
            FaceInterpolatedProperty<CoordinateTypes::Z>;

        struct HeatVolumetricCapacity
            : public StepPropertyGrid
        {
            static_assert(
                std::is_same<
                    Porosity::Grid_t,
                    SolidVolumetricHeatCapacity::Grid_t>::value
            );
            HeatVolumetricCapacity(
                const Porosity &porosity,
                const SolidVolumetricHeatCapacity &matrix_vol_heat_capacity,
                const Water &water)
                : StepPropertyGrid{
                    porosity.property_vals * water.volumetric_heat_capacity + (1 - porosity.property_vals) *matrix_vol_heat_capacity.property_vals,
                    porosity.grid
                }
            {}
        };

        template <typename Property_t, typename Grid_t>
        auto generate_log(
            const std::vector<RealType> &adata,
            const Grid_t &grid)
        {
            assert(grid.dual_stencils.size() == adata.size());

            return Property_t{
                Logs::StepProperty{
                    adata},
                grid};
        }

        struct HeatConductivity
            : public ZInterpolator,
              private AssertNonNegative
        {
            using ZInterpolator::FaceInterpolatedProperty;
        };

        // template <typename Property_t, typename Grid_t>
        // auto generate_faced_log(
        //     const std::vector<RealType> &adata,
        //     const Grid_t &grid)
        // {
        //     assert(grid.dual_stencils.size() == adata.size());

        //     return Property_t{
        //         Logs::StepProperty{
        //             adata},
        //         grid};
        // }

    } // Logs
} // GPN