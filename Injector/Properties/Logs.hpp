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
        /// It performs copying from the standard continer (Eigen or STL)
        /// to the local field variable.
        struct StepProperty
        {
            StepPropertyContainer data;

            StepProperty(
                const StepPropertyContainer &adata)
                : StepProperty{
                      std::vector<RealType>(
                          adata.cbegin(),
                          adata.cend())}
            {
            }

            StepProperty(
                const std::vector<RealType> &adata)
                : data(adata.size())
            {
#pragma region ASSERTIONS
                assert(data.size() > (decltype(data.size()))0);
                // for (auto idx{adata.cbegin()}; idx != adata.cend(); ++idx)
                //     // all properties are non-negative
                //     assert((*idx >= 0.0) || std::isnan(*idx));
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
            return StepProperty{StepPropertyContainer{v - rhs.data}};
        }

        struct StepPropertyGrid
        {
            using Grid_t = Grids::GridDual; // AxesGrid<CoordinateTypes::Z>;

            StepPropertyGrid(
                const StepPropertyContainer &interpolated_vals, // interpolated values corresponding to the grid
                const Grid_t &grid)
                : log_vals{interpolated_vals},
                  grid{grid}
            {
                assert(log_vals.size() == grid.dual_size() - 1ll);
                assert(log_vals.size() == grid.mesh_size());
                assert(grid.dual_stencils.dual_nodes.size() <= grid.dual_size());
            }

            StepPropertyGrid(const StepPropertyGrid &) noexcept = default;

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
                assert(property_vals.size() == grid.dual_stencils.size() - 1ll);
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

            operator const StepPropertyContainer &() const
            {
                return log_vals;
            }

            const StepPropertyContainer log_vals;
            const Grid_t grid;

        private:
            static StepPropertyContainer interpolate(
                const StepProperty &property_vals,
                const Grid_t &grid)
            {

                //    std::cout << "New interpolation:\n\n";

                StepPropertyContainer out(grid.mesh_size());
                // interpolate property_vals on the grid
                for (
                    // loop through all control_volumes
                    auto volume_id{0ll}, mesh_node_id{0ll};
                    volume_id < grid.dual_steps.size();
                    ++volume_id)
                {
                    //        std::cout << "volume_id =     " << volume_id << std::endl;
                    //        std::cout << "mesh_node_ids = " << std::endl;
                    // set constant value within a fixed control volume
                    for (;
                         (mesh_node_id < grid.dual_steps.size()) &&
                         (grid.mesh_nodes(mesh_node_id) < grid.dual_stencils(volume_id + 1ull));
                         ++mesh_node_id)
                    {
                        //                std::cout << mesh_node_id << ' ';
                        out(mesh_node_id) = property_vals.data(volume_id);
                    }
                    //        std::cout << std::endl;
                }

                return out;
            }
        };

        auto operator-(const StepPropertyGrid &lhs, RealType c)
        {
            return StepPropertyGrid{StepPropertyContainer{lhs.log_vals - c}, lhs.grid};
        }
        auto operator*(const StepPropertyGrid &lhs, RealType c)
        {
            return StepPropertyGrid{StepPropertyContainer{lhs.log_vals * c}, lhs.grid};
        }
        auto operator/(const StepPropertyGrid &lhs, RealType c)
        {
            return lhs * (1 / c);
        }
        auto operator*(const StepPropertyGrid &lhs, const StepPropertyGrid &rhs)
        {
            return StepPropertyGrid{StepPropertyContainer{lhs.log_vals * rhs.log_vals}, lhs.grid};
        }
        auto operator/(const StepPropertyGrid &lhs, const StepPropertyGrid &rhs)
        {
            return StepPropertyGrid{StepPropertyContainer{lhs.log_vals / rhs.log_vals}, lhs.grid};
        }

        auto operator+(const StepPropertyGrid &lhs, const StepPropertyGrid &rhs)
        {
            return StepPropertyGrid{StepPropertyContainer{(lhs.log_vals + rhs.log_vals)}, lhs.grid};
        }

        auto operator-(RealType c, const StepPropertyGrid &rhs)
        {
            return StepPropertyGrid{StepPropertyContainer{c - rhs.log_vals}, rhs.grid};
        }

        struct AssertNonNegative
        {
            AssertNonNegative(const StepPropertyGrid &vals)
            {
                const auto &data{vals.log_vals};
                assert(std::all_of(
                    data.cbegin(), data.cend(),
                    [](const RealType x)
                    { return x >= 0.0; }));
            }
        };
#pragma region INDICATORS
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

        /// @brief Indicator of permeable layers,
        /// so the liquid flows along these layers
        struct IsPermeable : public IndicatorProperty
        {
            using IndicatorProperty::IndicatorProperty;
        };

        /// @brief Indicator of perforated cells,
        /// so the liquid can leave the column,
        /// to further flow across the cement to the permeable layers
        struct IsPerforated : public IndicatorProperty
        {
            using IndicatorProperty::IndicatorProperty;
        };

        /// @brief Indicator of permeable layers
        /// in the reservoir,
        /// which can accept fluid
        struct IsGhostLayer : public IndicatorProperty
        {
            using IndicatorProperty::IndicatorProperty;
        };

        /// @brief Indicator of cells with damaged column,
        /// so the liquid can leave the column,
        /// to further flow vertically along the cement.
        /// Ignores perforated cells.
        struct IsDamaged : public IndicatorProperty
        {
            using IndicatorProperty::IndicatorProperty;
        };
#pragma endregion

#pragma region HYDRODYNAMICS-LOGS
        /// @brief Pressure at the external boundary of the
        /// computation domain. far from the well
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
                // for (std::ptrdiff_t id{0ll}; id < pressure.size(); ++id)
                //     assert(
                //         ((is_permeable(id) == 1.0) && (pressure(id) > 0.0)) ||
                //         ((is_permeable(id) == 0.0) && (pressure(id) == 0.0)));
            }
        };

        struct HydrostaticPressure
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            HydrostaticPressure(
                const StepPropertyGrid &pressure)
                : StepPropertyGrid{pressure},
                  AssertNonNegative{pressure}
            {
            }
        };

        namespace InternalUse
        {
            /// @brief Rate distribution along the
            /// layers
            struct RateWeights
                : public StepPropertyGrid,
                  private AssertNonNegative
            {
                RateWeights(
                    const StepPropertyGrid &weights,
                    const StepPropertyGrid &indicator)
                    : StepPropertyGrid{normalize(weights)},
                      AssertNonNegative{weights}
                {
                    assert(weights.size() == indicator.size());
                    for (std::ptrdiff_t id{0ll}; id < weights.size(); ++id)
                        assert(
                            ((indicator(id) == 1.0)) ||
                            ((indicator(id) == 0.0) && (weights(id) == 0.0)));

                    assert(
                        std::all_of(weights.log_vals.cbegin(), weights.log_vals.cend(), [](const RealType v)
                                    { return v == 0.0; }) ||
                        (std::abs(log_vals.sum() - 1.0) < 1E-12));
                }

            private:
                static StepPropertyGrid normalize(const StepPropertyGrid &weights)
                {
                    const RealType sum{weights.log_vals.sum()};
                    return StepPropertyGrid{weights.log_vals / (sum == 0.0 ? 1.0 : sum), weights.grid};
                }
            };
        } // InternalUse

        struct RFP
            : public StepPropertyGrid
        {
            RFP(const StepPropertyGrid &rfp,
                        const IsPermeable &indicator)
                        :StepPropertyGrid{rfp*indicator}
            {
                assert(rfp.size() == indicator.size());
                for (std::ptrdiff_t id{0ll}; id < rfp.size(); ++id)
                    assert(
                        ((indicator(id) == 1.0)) ||
                        ((indicator(id) == 0.0) && (rfp(id) == 0.0)));
            }
        };
        
        struct RFP_weights
            : public RFP
        {
            RFP_weights(const StepPropertyGrid &rfp,
                        const IsPermeable &is_permeable)
                : RFP{ InternalUse::RateWeights(rfp, is_permeable), is_permeable}
            {
            }
        };

        struct WFP
            : public StepPropertyGrid
        {
            WFP(const StepPropertyGrid &wfp,
                        const IsPerforated &indicator)
                : StepPropertyGrid{wfp*indicator}
            {
                assert(wfp.size() == indicator.size());
                for (std::ptrdiff_t id{0ll}; id < wfp.size(); ++id)
                    assert(
                        ((indicator(id) == 1.0)) ||
                        ((indicator(id) == 0.0) && (wfp(id) == 0.0)));
            }
        };
        
        struct WFP_weights
            : public WFP
        {
            WFP_weights(const StepPropertyGrid &wfp,
                        const IsPerforated &is_perforated)
                : WFP{InternalUse::RateWeights{wfp, is_perforated}, is_perforated}
            {
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

        struct MediumCompressibility
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            MediumCompressibility(
                const StepPropertyGrid &compressibility,
                const IsPermeable &is_permeable)
                : StepPropertyGrid{compressibility},
                  AssertNonNegative{compressibility}
            {
                assert(compressibility.size() == is_permeable.size());
                // for (std::ptrdiff_t id{0ll}; id < compressibility.size(); ++id)
                //     assert(
                //         ((is_permeable(id) == 1.0) && (compressibility(id) > 0.0)) ||
                //         ((is_permeable(id) == 0.0) && (compressibility(id) == 0.0)));
            }
        };

        struct SkinFactor
            : public StepPropertyGrid
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
#pragma endregion
#pragma region HEAT-LOGS
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

        struct Geotherma
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            Geotherma(
                const StepPropertyGrid &temperature)
                : StepPropertyGrid{temperature},
                  AssertNonNegative{temperature}
            {
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

        struct MediumHeatVolumetricCapacity
            : public StepPropertyGrid
        {
            MediumHeatVolumetricCapacity(
                const StepPropertyGrid &capacity)
                : StepPropertyGrid{capacity}
            {
            }
        };

        struct HeatConductivity
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            HeatConductivity(
                const StepPropertyGrid &conductivity)
                : StepPropertyGrid{conductivity},
                  AssertNonNegative{conductivity}
            {
            }
        };

        using MediumHeatConductivity = HeatConductivity;

        struct ThermalDiffusivity
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            ThermalDiffusivity(
                const StepPropertyGrid &conductivity)
                : StepPropertyGrid{conductivity},
                  AssertNonNegative{conductivity}
            {
            }
        };
#pragma endregion
        template <typename Property_t, typename Grid_t>
        [[deprecated]]
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

    } // Logs
} // GPN