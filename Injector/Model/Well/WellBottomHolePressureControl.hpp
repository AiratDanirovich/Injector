#pragma once

#include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/Model/Well/DefaultWellNmerics.hpp>
#include <Injector/Model/Well/SomeWell.hpp>

#include <Injector/Model/Collector.hpp>
#include <Injector/Properties/Logs.hpp>

#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>

namespace GPN
{
    namespace Wells
    {
        namespace BotHolePresControl
        {
#pragma region BOUNDARY-CONDITION
            template <typename Fluid_t, typename History_t, typename Grid2D_t>
            struct BotHolePresFunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
            {
            protected:
                using HydrostaticPressureFactory_t = Logs::HydrostaticPressureFactory<Fluid_t, typename Grid2D_t::Axes1Coordinate_t>;

            public:
                /// @brief
                /// @param history
                /// @param ext_pressure
                /// @param PI productivity index that multiplies pressure difference at sandface
                /// @param grid_ptr
                BotHolePresFunctorBC(
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

                        const auto sandface_pressure{hydrostatic_factory.create(history->pressure()).log_vals};
                        const auto out{sandface_pressure(z_id)};
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

                const HydrostaticPressureFactory_t &hydrostatic_factory;
            };

            struct BotHolePresBC : public BoundaryConditions::GeneralBC
            {
                template <typename Grid2D_t>
                BotHolePresBC(const cptr<Grid2D_t> &grid,
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
            struct WellBottomHolePressureControl
                : public DefaultWellNumerics<0ll>,
                  public SomeWell<History_t, Grid2D_t, CrossFlow_t>
            {
                using Base = SomeWell<History_t, Grid2D_t, CrossFlow_t>;

                using functor_type = BotHolePresFunctorBC<Fluid_t, History_t, Grid2D_t>;
                using hydro_bc_type = BotHolePresBC;
                using grid_type = Grid2D_t;

                template <
                    typename Grid2D_t,
                    typename Capacity_t,
                    typename ConvectionTermFactory_t,
                    typename BC_t>
                using solver_type =
                    EqSolver::FullImplicit::Solver<
                        Grid2D_t, Capacity_t, ConvectionTermFactory_t, BC_t>;

                const ptr<const hydro_bc_type> hydro_bc;

                WellBottomHolePressureControl(
                    const Properties::Rocks::RocksProps<Grid2D_t> &
                        rock_field_props,
                    const CrossFlow_t &well_base,
                    const cptr<History_t> history,
                    const Fluid_t &fluid,
                    const cptr<Grid2D_t> grid2D_rocks)
                    : Base{rock_field_props, well_base, history, grid2D_rocks},
                      size{grid2D_rocks->first_coord().mesh_size()},
                      fluid{fluid},
                      hydrostatic_factory{history->z_ref, fluid, grid2D_rocks->first_coord()}
                {
                    const_cast<ptr<const hydro_bc_type> &>(hydro_bc) =
                        std::make_shared<const hydro_bc_type>(
                            grid2D_rocks,
                            std::make_shared<const functor_type>(
                                history, rock_field_props.base_hydrodynamics.ext_pressure,
                                hydrostatic_factory,
                                Base::PI,
                                grid2D_rocks));
                }

                template <typename HistoryRecord_t>
                RealType get_total_bottomhole_rate(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    return get_RFP(record).log_vals.sum();
                }

                template <typename HistoryRecord_t>
                void set_well_flow_field(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure)
                {
                    const auto well_pres_prof{
                        well_pressure_profile(record, collector_pressure).log_vals};
                    Base::set_flux(
                        Base::set_RFP(
                            well_pres_prof, collector_pressure));
                    Base::set_well_flow_field(
                        Base::get_verticle_well_flow(record),
                        Base::get_verticle_cement_flow(record),
                        Base::get_WFP(record), Base::get_RFP(record));

                                        const StepPropertyContainer depression_at_sandface{
                        -(collector_pressure.col(0ll) - record.pressure * Base::is_permeable.log_vals)};
                    assert(
                        std::all_of(
                            depression_at_sandface.cbegin(),
                            depression_at_sandface.cend(),
                            [](const RealType v)
                            {
                                return !std::isnan(v);
                            }));
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
                    return well_pressure_profile(record, collector_pressure).log_vals * Base::is_permeable.log_vals;
                }

                const Fluid_t &fluid;
                const std::ptrdiff_t size;
                const Logs::HydrostaticPressureFactory<Fluid_t, typename Grid2D_t::Axes1Coordinate_t>
                    hydrostatic_factory;

                template <typename HistoryRecord_t>
                const RealType P_bot(
                    const HistoryRecord_t &record,
                    const auto & /*collector_pressure*/) const
                {
                    return record.pressure;
                }
            };
        } // BotHolePresControl
    } // Wells
} // GPN