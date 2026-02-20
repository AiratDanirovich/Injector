#pragma once

#include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/Model/Well/DefaultWellNmerics.hpp>
#include <Injector/Model/Well/SomeWell.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Properties/Logs.hpp>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/FullImplicit/WellRateControlSolver.hpp>

namespace GPN
{
    namespace Wells
    {
        namespace BotHoleRateControl
        {
#pragma region BOUNDARY-CONDITION
            template <typename Fluid_t, typename History_t, typename Grid2D_t>
            struct BotHoleRateFunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
            {
            protected:
                using HydrostaticPressureFactory_t = Logs::HydrostaticPressureFactory<Fluid_t, typename Grid2D_t::Axes1Coordinate_t>;

            public:
                /// @brief
                /// @param history
                /// @param ext_pressure
                /// @param PI productivity index that multiplies pressure difference at sandface
                /// @param grid_ptr
                BotHoleRateFunctorBC(
                    //    const RealType fluid_density,
                    const ptr<const History_t> history,
                    const Logs::ExternalPressure &ext_pressure,
                    const HydrostaticPressureFactory_t &hydrostatic_factory,
                    const StepPropertyContainer &PI,
                    const ptr<const Grid2D_t> grid_ptr)
                    : BoundaryConditions::GeneralBC::BCFunctorBase{},
                      history{history},
                      ext_pressure{ext_pressure},
                      PI{PI},
                      grid_ptr{grid_ptr},
                      fluid_weight{hydrostatic_factory.fluid.density * Gravity::value()},
                      hydrostatic_factory{hydrostatic_factory}
                {
                    static_assert(Grid2D_t::l_margin == 3ll);
                }

                BC_descriptor operator()(const ptrdiff_t z_id, const RealType r, const RealType,
                                         const BCType bc_type) const override
                {
                    // boundary condition at sandface
                    if (r == grid_ptr->second_coord().dual_front())
                    {
                        assert(bc_type == BCType::third);
                        // const auto out{history->pressure()};
                        const auto dz{hydrostatic_factory.grid_z.mesh_nodes(z_id) - hydrostatic_factory.z_ref};
                        const auto out{fluid_weight * dz};
                        return BC_descriptor::BC_III(out, PI(z_id));
                        //    return out;
                    }

                    // boundary condition at external contour
                    if (r == grid_ptr->second_coord().dual_back())
                    {
                        assert(bc_type == BCType::first);
                        return BC_descriptor::BC_I(ext_pressure(z_id));
                        //    return ext_pressure(z_id);
                    }

                    assert(false);
                    throw std::runtime_error("Wrong radial coordinate in BC.");
                    //    return 0.0;
                }

                BC_descriptor operator()(const RealType z, const ptrdiff_t, const RealType,
                                         const BCType bc_type) const override
                {
                    // boundary conditions are set exactly at domain boundaries
                    assert((z == grid_ptr->first_coord().dual_front()) || (z == grid_ptr->first_coord().dual_back()));
                    // assume zero "diffusion" flux in hydrodynamic equation
                    assert(bc_type == BCType::second);
                    return BC_descriptor::BC_II(0.0);
                    //    return 0.0;
                }

                const StepPropertyContainer &PI;
                const ptr<const History_t> history;

            protected:
                const Logs::ExternalPressure &ext_pressure;
                const ptr<const Grid2D_t> grid_ptr;
                const RealType fluid_weight;

                const HydrostaticPressureFactory_t &hydrostatic_factory;
            };

            struct BotHoleRateBC : public BoundaryConditions::GeneralBC
            {
                template <typename Grid2D_t>
                BotHoleRateBC(const cptr<Grid2D_t> &grid,
                              const cptr<const BCFunctorBase> functor)
                    : BoundaryConditions::GeneralBC{grid, functor, BoundaryCondition::BCType::second}
                {
                    bc_types[west_id] = BoundaryCondition::BCType::third;
                    bc_types[east_id] = BoundaryCondition::BCType::first;
                    east_bc_type = std::vector<BoundaryCondition::BCType>(grid->first_coord().mesh_size(), bc_types[east_id]);
                }
                void set_bc_type(const RealType t)
                {
                    this->t = t;
                }
            };
#pragma endregion
            template <
                typename History_t,
                typename Fluid_t,
                typename Grid2D_t,
                typename CrossFlow_t>
            struct WellBottomHoleRateControl
                : public DefaultWellNumerics<1ll>,
                  public SomeWell<History_t, Grid2D_t, CrossFlow_t>
            {
                using Base = SomeWell<History_t, Grid2D_t, CrossFlow_t>;

                using functor_type = BotHoleRateFunctorBC<Fluid_t, History_t, Grid2D_t>;
                using hydro_bc_type = BotHoleRateBC;
                using grid_type = Grid2D_t;

                template <
                    typename Grid2D_t,
                    typename Capacity_t,
                    typename ConvectionTermFactory_t,
                    typename BC_t>
                using solver_type =
                    EqSolver::FullImplicit::WellRateControlSolver<
                        Grid2D_t, Capacity_t, ConvectionTermFactory_t, BC_t,
                        WellBottomHoleRateControl<History_t, Fluid_t, Grid2D_t, CrossFlow_t>>;

                const ptr<const hydro_bc_type> hydro_bc;

                WellBottomHoleRateControl(
                    const Properties::Rocks::RocksProps<Grid2D_t> &
                        rock_field_props,
                    const CrossFlow_t &well_base,
                    const cptr<History_t> history,
                    const Fluid_t &fluid,
                    const cptr<Grid2D_t> grid2D_rocks)
                    : Base{rock_field_props, well_base, history, grid2D_rocks},
                      size{grid2D_rocks->first_coord().mesh_size()},
                      is_permeable{rock_field_props.base_hydrodynamics.is_permeable},
                      fluid{fluid},
                      hydrostatic_factory{history->z_ref, fluid, grid2D_rocks->first_coord()}
                {
                    const_cast<ptr<const hydro_bc_type> &>(hydro_bc) =
                        std::make_shared<const hydro_bc_type>(
                            grid2D_rocks,
                            std::make_shared<const functor_type>(
                                history,
                                rock_field_props.base_hydrodynamics.ext_pressure,
                                hydrostatic_factory,
                                Base::PI,
                                grid2D_rocks));
                }

                template <typename HistoryRecord_t>
                RealType get_total_bottomhole_rate(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    return record.rate;
                }

                template <typename HistoryRecord_t>
                void set_well_flow_field(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure)
                {
                    Base::set_flux(set_RFP(record, collector_pressure));
                    Base::set_well_flow_field(
                        Base::get_verticle_well_flow(record),
                        Base::get_verticle_cement_flow(record),
                        Base::get_WFP(record), Base::get_RFP(record));
                }

                template <typename HistoryRecord_t>
                const auto well_pressure_profile(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    return hydrostatic_factory.create(
                        P_bot(Base::history->get_current_record(), collector_pressure));
                }

                template <typename HistoryRecord_t>
                StepPropertyContainer get_pressure_at_sandface(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    // Pressure at the level of NON-permeable layers is assumed zero.
                    // Pressure at permeable layers is equal to history->pressure()
                    return well_pressure_profile(record, collector_pressure).log_vals * is_permeable.log_vals;
                }

                const Fluid_t &fluid;
                const std::ptrdiff_t size;
                const Logs::IsPermeable &is_permeable;
                const Logs::HydrostaticPressureFactory<Fluid_t, typename Grid2D_t::Axes1Coordinate_t>
                    hydrostatic_factory;

                template <typename HistoryRecord_t>
                const RealType P_bot(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    const auto rate{record.rate};
                    const auto weight{fluid.density * Gravity::value()};
                    const auto dz{Base::grid2D_rocks->first_coord().mesh_nodes - Base::z_ref()};
                    const RealType term1{(Base::PI * (weight * dz - collector_pressure.col(0ll))).sum()};
                    assert(!std::isnan(rate) && !std::isinf(rate));
                    const auto out{(rate - term1) / Base::PI.sum()};
                    return out;
                }

            protected:
                template <typename HistoryRecord_t>
                const auto set_RFP(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    return Logs::RFPFactory::create_from_container<Logs::RFP>(
                        StepPropertyContainer{(
                                                  Base::PI * (well_pressure_profile(record, collector_pressure).log_vals -
                                                              collector_pressure.col(0ll)))
                                                  .eval()},
                        is_permeable);
                }
            };
        } // BotHoleRateControl
    } // Wells
} // GPN