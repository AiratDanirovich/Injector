#pragma once

#include <algorithm>

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

        struct InterpolatedDataContainer : private StepPropertyContainer
        {
            using StepPropertyContainer::StepPropertyContainer;
            using StepPropertyContainer::operator*;
            using StepPropertyContainer::operator/;
            using StepPropertyContainer::operator+;
            using StepPropertyContainer::operator();
            using StepPropertyContainer::size;
        };

        /// @brief Container for values of step properties.
        /// Copies from standard continer (Eigen or STL)
        /// to local field variable.
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
            {
#pragma region ASSERTIONS
                assert(data.size() > (decltype(data.size()))0);
                for (auto idx{adata.cbegin()}; idx != adata.cend(); ++idx)
                    // all properties are non-negative
                    assert(*idx >= 0.0);
#pragma endregion
                std::copy(adata.cbegin(), adata.cend(), data.begin());
            }

            StepProperty(const StepProperty &) noexcept = default;

            auto size() const
            {
                return data.size();
            }

            auto operator()(auto id) const
            {
                return data(id);
            }
        };

        auto operator-(RealType v, const StepProperty &rhs)
        {
            return StepProperty{v - rhs.data};
        }

        struct StepPropertyGrid
        {
            using Grid_t = Grids::GridDual; // AxesGrid<CoordinateTypes::Z>;

            StepPropertyGrid(const StepPropertyGrid &) noexcept = default;
            StepPropertyGrid(
                const StepPropertyContainer &interpolated_vals, // interpolated values corresponding to the grid
                const Grid_t &grid)
                : log_vals{interpolated_vals},
                  grid{grid}
            {
                assert(log_vals.size() == grid.dual_size() - 1ll);
            }

            StepPropertyGrid(
                // stencil values per every layer
                const StepProperty &property_vals,
                const Grid_t &grid)
                : StepPropertyGrid{
                      // interpolated stencil-values
                      // over refined grid nodes
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
                return grid.dual_stencils.size();
            }

            auto size() const
            {
                return log_vals.size();
            }

            auto operator-(RealType c) const
            {
                const auto &lhs{*this};
                return StepPropertyGrid{lhs.log_vals - c, lhs.grid};
            }
            auto operator*(RealType c) const
            {
                const auto &lhs{*this};
                return StepPropertyGrid{lhs.log_vals * c, lhs.grid};
            }
            auto operator*(const StepPropertyGrid &rhs) const
            {
                const auto &lhs{*this};
                return StepPropertyGrid{lhs.log_vals * rhs.log_vals, lhs.grid};
            }
            auto operator/(const StepPropertyGrid &rhs) const
            {
                const auto &lhs{*this};
                return StepPropertyGrid{lhs.log_vals / rhs.log_vals, lhs.grid};
            }

            auto operator+(const StepPropertyGrid &rhs) const
            {
                const auto &lhs{*this};
                return StepPropertyGrid{(lhs.log_vals + rhs.log_vals), lhs.grid};
            }

            const StepPropertyContainer log_vals;
            const Grid_t grid;

        private:
            static StepPropertyContainer interpolate(
                const StepProperty &property_vals,
                const Grid_t &grid)
            {
                StepPropertyContainer out(grid.mesh_size());
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

        struct AssertNonNegative
        {
            AssertNonNegative(const StepPropertyGrid &vals)
            {
                const auto &data{vals.log_vals};
                std::for_each(data.cbegin(), data.cend(),
                              [](RealType x)
                              { assert(x >= 0.0); });
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
                const auto &data{vals.log_vals};
                std::for_each(data.cbegin(), data.cend(),
                              [](RealType x)
                              { assert(x == 0.0 || x == 1.0); });
            }
        };

        struct IsPermeable : public IndicatorProperty
        {
            using IndicatorProperty::IndicatorProperty;
        };

        struct ExternalPressure
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            ExternalPressure(
                const StepPropertyGrid &pressure,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{pressure},
                  AssertNonNegative{pressure}
            {
                assert(pressure.size() == is_permeable.size());
                for (std::ptrdiff_t id{0ll}; id < pressure.size(); ++id)
                    assert(
                        ((is_permeable(id) == 1.0) && (pressure(id) > 0.0)) ||
                        ((is_permeable(id) == 0.0) && (pressure(id) == 0.0)));
            }
        };

        struct RFP
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            RFP(const StepPropertyGrid &rfp,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{rfp},
                  AssertNonNegative{rfp}
            {
                assert(rfp.size() == is_permeable.size());
                for (std::ptrdiff_t id{0ll}; id < rfp.size(); ++id)
                    assert(
                        ((is_permeable(id) == 1.0) && (rfp(id) > 0.0)) ||
                        ((is_permeable(id) == 0.0) && (rfp(id) == 0.0)));
            }
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
                const StepPropertyGrid &density)
                : StepPropertyGrid{density},
                  AssertNonNegative{density}
            {
            }
        };
        struct SolidSpecificHeatCapacity
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            SolidSpecificHeatCapacity(
                const StepPropertyGrid &capacity)
                : StepPropertyGrid{capacity},
                  AssertNonNegative{capacity}
            {
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
                      density * heat_capacity}
            {
            }
        };

        struct HeatVolumetricCapacity
            : public StepPropertyGrid
        {
            static_assert(
                std::is_same<
                    Porosity::Grid_t,
                    SolidVolumetricHeatCapacity::Grid_t>::value);
            HeatVolumetricCapacity(
                const Porosity &porosity,
                const SolidVolumetricHeatCapacity &matrix_vol_heat_capacity,
                const Water &water)
                : StepPropertyGrid{
                      porosity * water.volumetric_heat_capacity + (porosity - 1.0) * (-1.0) * matrix_vol_heat_capacity}
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
                // if dual_size() == 2 --- no internal faces, only single cell
                assert(grid.dual_size() > 2ll);

                InternalFaceValues out(log.grid.dual_size() - 2ll);

                for (auto id{0ll}; id < out.size(); ++id)
                    out(id) = Axes::face_interpolator(
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
            HeatConductivity(
                const StepPropertyGrid &conductivity)
                : ZInterpolator{conductivity},
                  AssertNonNegative{conductivity}
            {
            }
        };

        struct ThermalDiffusivity
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            ThermalDiffusivity(
                const HeatVolumetricCapacity &capacity,
                const HeatConductivity &conductivity)
                : StepPropertyGrid{
                      conductivity / capacity},
                  AssertNonNegative{capacity}
            {
            }
        };

    } // Logs
} // GPN