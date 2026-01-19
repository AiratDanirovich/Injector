#pragma once

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/BoundaryConditions.hpp>

namespace GPN
{
    namespace Wells
    {
        namespace BotHolePresControl
        {
#pragma region BOUNDARY-CONDITION
            template <typename History_t, typename Grid2D_t>
            struct FunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
            {
                FunctorBC(
                    const cptr<const History_t> history,
                    const Logs::ExternalPressure &ext_pressure,
                    const cptr<const Grid2D_t> grid_ptr)
                    : history{history},
                      ext_pressure{ext_pressure},
                      grid_ptr{grid_ptr}
                {
                }

                RealType operator()(const ptrdiff_t z_id, const RealType r, const RealType,
                                    const BCType bc_type) const override
                {
                    // boundary condition at sandface
                    if (r == grid_ptr->second_coord().dual_front())
                    {
                        assert(bc_type == BCType::first);
                        const auto out{history->pressure()};
                        return out;
                    }

                    // boundary condition at external contour
                    if (r == grid_ptr->second_coord().dual_back())
                    {
                        assert(bc_type == BCType::first);
                        return ext_pressure(z_id);
                    }

                    assert(false);
                    return 0.0;
                }

                RealType operator()(const RealType z, const ptrdiff_t, const RealType,
                                    const BCType bc_type) const override
                {
                    // boundary conditions are set exactly at domain boundaries
                    assert((z == grid_ptr->first_coord().dual_front()) || (z == grid_ptr->first_coord().dual_back()));
                    // assume zero "diffusion" flux in hydrodynamic equation
                    assert(bc_type == BCType::second);
                    return 0.0;
                }

            protected:
                const Logs::ExternalPressure &ext_pressure;
                const cptr<const Grid2D_t> grid_ptr;
                const cptr<const History_t> history;
            };

            struct BotHolePresBC : public BoundaryConditions::GeneralBC
            {
                template <typename Grid2D_t>
                BotHolePresBC(const cptr<Grid2D_t> &grid,
                              const cptr<const BCFunctorBase> functor)
                    : BoundaryConditions::GeneralBC{grid, functor, BoundaryCondition::BCType::first}
                {
                    bc_types[north_id] = BoundaryCondition::BCType::second;
                    bc_types[south_id] = BoundaryCondition::BCType::second;
                }
                void set_bc_type(const RealType t)
                {
                    this->t = t;
                }
            };
#pragma endregion

            template <
                typename History_t,
                typename Grid2D_t,
                typename Well_t>
            struct WellBottomHolePressureControl : public Well_t
            {
                using functor_type = FunctorBC<History_t, Grid2D_t>;

                using hydro_bc_type = BotHolePresBC;

                const ptr<const hydro_bc_type> hydro_bc;

                WellBottomHolePressureControl(
                    const Properties::Rocks::RocksProps<Grid2D_t> &
                        rock_field_props,
                    const Well_t &well_base,
                    const cptr<History_t> history,
                    const cptr<Grid2D_t> grid2D_rocks)
                    : Well_t{well_base},
                      hydro_bc{
                          std::make_shared<const hydro_bc_type>(
                              grid2D_rocks,
                              std::make_shared<const functor_type>(
                                  history, rock_field_props.base_hydrodynamics.ext_pressure,
                                  grid2D_rocks))},
                      size{grid2D_rocks->first_coord().mesh_size()},
                      mobility{set_mobility(rock_field_props, grid2D_rocks)}
                {
                }

                template <typename HistoryRecord_t>
                StepPropertyContainer get_pressure_at_symmetry_axis(
                    const HistoryRecord_t &record,
                    const auto &ref_pressure) const
                {
                    return StepPropertyContainer::Constant(size, record.pressure);
                }

                // template <typename HistoryRecord_t>
                // StepPropertyContainer get_RFP(
                //     const HistoryRecord_t &record,
                //     const auto &ref_pressure) const
                // {
                //     return mobility * (ref_pressure.col(0ll) - record.pressure);
                // }
                
                template <typename HistoryRecord_t>
                RealType get_total_bottomhole_rate(
                    const HistoryRecord_t &record,
                    const auto &ref_pressure) const
                {
                    return get_RFP(record, ref_pressure).sum();
                }


                template <typename HistoryRecord_t>
                StepPropertyContainer get_total_bottomhole_rate(
                    const HistoryRecord_t &record,
                    const auto &ref_pressure) const
                {
                    return mobility * (ref_pressure.col(0ll) - record.pressure);
                }

                const std::ptrdiff_t size;

            protected:
                static auto set_mobility(
                    const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                    cptr<Grid2D_t> grid2D_rocks)
                {
                    // center of cell next to the well
                    const auto r3{grid2D_rocks->second_coord().mesh_nodes(0ll)};
                    // sandface
                    const auto r2_face{grid2D_rocks->second_coord().dual_nodes(0ll)};
                    const auto two_pi{2.0 * std::numbers::pi_v<RealType>};
                    const auto is_permeable{rock_field_props.base_hydrodynamics.is_permeable.log_vals};

                    Eigen::ArrayX<RealType> out{two_pi / std::log(r3 / r2_face) *
                                                rock_field_props.mobility_axes2.col(Grid2D_t::l_margin) *
                                                grid2D_rocks->first_coord().volumes()};

                    for (auto row{0ll}; row < grid2D_rocks->first_coord().mesh_size(); ++row)
                    {
                        if (is_permeable(row) == 0.0)
                            out(row) = 0.0;
                    }
                    return out;
                }

                const StepPropertyContainer mobility;
            };

        } // BotHolePresControl
    } // Wells
} // GPN