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
        namespace ResFlowProfileControl
        {
#pragma region BOUNDARY-CONDITION
            template <typename History_t, typename Grid2D_t>
            struct RFPControlFunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
            {
                /// @brief
                /// @param history
                /// @param ext_pressure
                /// @param rfp log_vals of Logs::RFP_weights container
                /// @param grid_ptr Grid in collector, outside the sandface
                RFPControlFunctorBC(
                    const ptr<const History_t> history,
                    const Logs::ExternalPressure &ext_pressure,
                    const StepPropertyContainer &rfp,
                    const Logs::IsPermeable &is_permeable,
                    const ptr<const Grid2D_t> grid_ptr)
                    : BoundaryConditions::GeneralBC::BCFunctorBase{},
                      history{history},
                      ext_pressure{ext_pressure},
                      grid_ptr{grid_ptr},
                      rfp{Logs::RFPFactory::create_from_container<Logs::RFP_weights>(
                          rfp, is_permeable)}
                {
                    static_assert(Grid2D_t::l_margin == 3ll);
                    assert(std::abs(this->rfp.log_vals.sum() - 1.0) < 1e-12);
                }

                BC_descriptor operator()(const ptrdiff_t z_id, const RealType r, const RealType,
                                         const BCType bc_type) const override
                {
                    // sandface
                    if (r == grid_ptr->second_coord().dual_front())
                    {
                        assert(bc_type == BCType::second);
                        const auto out{rfp(z_id) * history->rate()};
                        return BC_descriptor::BC_II(out);
                        //    return out;
                    }
                    // external contour
                    if (r == grid_ptr->second_coord().dual_back())
                    {
                        assert(bc_type == BCType::first);
                        return BC_descriptor::BC_I(ext_pressure(z_id));
                        //    return ext_pressure(z_id);
                    }

                    assert(false);
                    throw std::runtime_error("Wrong radial coordinate in BC.");
                }

                BC_descriptor operator()(const RealType z, const ptrdiff_t, const RealType,
                                         const BCType bc_type) const override
                {
                    // boundary conditions are set exactly at domain boundaries
                    assert((z == grid_ptr->first_coord().dual_front()) || (z == grid_ptr->first_coord().dual_back()));
                    // assume zero diffusion flux in hydrodynamic equation
                    assert(bc_type == BCType::second);
                    return BC_descriptor::BC_II(0.0);
                    //    return 0.0;
                }

            protected:
                const Logs::ExternalPressure &ext_pressure;
                const ptr<const Grid2D_t> grid_ptr;
                const Logs::RFP_weights rfp;
                const ptr<const History_t> history;
            };

            struct RFPControlBC : public BoundaryConditions::GeneralBC
            {
                template <typename Grid2D_t>
                RFPControlBC(const cptr<Grid2D_t> &grid,
                             const ptr<const BCFunctorBase> functor)
                    : BoundaryConditions::GeneralBC{grid, functor, BoundaryCondition::BCType::second}
                {
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
            struct WellReservoirFlowProfileControl
                : 
                  public DefaultWellNumerics<0ll>,
                  public SomeWell<Grid2D_t, CrossFlow_t>
            {
                using Base = SomeWell<Grid2D_t, CrossFlow_t>;

                using functor_type = RFPControlFunctorBC<History_t, Grid2D_t>;
                using hydro_bc_type = RFPControlBC;
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

                WellReservoirFlowProfileControl(
                    const Properties::Rocks::RocksProps<Grid2D_t> &
                        rock_field_props,
                    const CrossFlow_t &well_base,
                    const cptr<History_t> history,
                    const Fluid_t &,
                    const cptr<Grid2D_t> grid2D_rocks)
                    : 
                      Base{well_base, grid2D_rocks},
                      resistivity{set_resistivity(rock_field_props, well_base, grid2D_rocks)},
                      hydro_bc{
                          std::make_shared<const hydro_bc_type>(
                              grid2D_rocks,
                              std::make_shared<const functor_type>(
                                  history, rock_field_props.base_hydrodynamics.ext_pressure,
                                  well_base.rfp, rock_field_props.base_hydrodynamics.is_permeable,
                                  grid2D_rocks))},
                      history{history},
                      rock_field_props{rock_field_props}
                {
                    assert(std::abs(well_base.rfp.sum() - 1.0) < 1e-12);
                    assert(std::abs(well_base.wfp.sum() - 1.0) < 1e-12);
                }

                template <typename HistoryRecord_t>
                const auto get_RFP(const HistoryRecord_t &record) const
                {
                    return Base::rfp * record.rate;
                }
                template <typename HistoryRecord_t>
                const auto get_WFP(const HistoryRecord_t &record) const
                {
                    return Base::wfp * record.rate;
                }
                template <typename HistoryRecord_t>
                const auto get_verticle_cement_flow(const HistoryRecord_t &record) const
                {
                    return Base::verticle_flux_in_cement * record.rate;
                }
                template <typename HistoryRecord_t>
                const auto get_verticle_well_flow(const HistoryRecord_t &record) const
                {
                    return Base::verticle_flux_in_well * record.rate;
                }

                template <typename HistoryRecord_t>
                StepPropertyContainer get_pressure_at_sandface(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    return collector_pressure.col(0ll) + record.rate * resistivity;
                }

                template <typename HistoryRecord_t>
                RealType get_total_bottomhole_rate(
                    const HistoryRecord_t &record,
                    const auto & /*collector_pressure*/) const
                {
                    return record.rate;
                }

                template <typename HistoryRecord_t>
                void set_well_flow_field(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure)
                {
                    const auto rfp_{get_RFP(record)};
                    const auto wfp_{get_WFP(record)};
                    Base::set_well_flow_field(wfp_, rfp_);

                    // set verticle flux
                    Base::flow_axes1_value.col(0ll) = get_verticle_well_flow(record);
                    Base::flow_axes1_value.col(1ll) = 0.0;
                    Base::flow_axes1_value.col(2ll) = get_verticle_cement_flow(record);
                }

                const cptr<History_t> history;
                const Properties::Rocks::RocksProps<Grid2D_t> &
                    rock_field_props;

            private:
                static auto set_resistivity(
                    const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                    const CrossFlow_t &well,
                    cptr<Grid2D_t> grid2D_rocks)
                {
                    const auto r3{grid2D_rocks->second_coord().mesh_nodes(0ll)};
                    const auto r2_face{grid2D_rocks->second_coord().dual_nodes(0ll)};
                    const auto two_pi{2.0 * std::numbers::pi_v<RealType>};
                    const auto is_permeable{rock_field_props.base_hydrodynamics.is_permeable.log_vals};

                    Eigen::ArrayX<RealType> out{well.rfp / (two_pi / std::log(r3 / r2_face) *
                                                            rock_field_props.mobility_axes2.col(Grid2D_t::l_margin) *
                                                            grid2D_rocks->first_coord().volumes())};

                    for (auto row{0ll}; row < grid2D_rocks->first_coord().mesh_size(); ++row)
                    {
                        if (is_permeable(row) == 0.0)
                            out(row) = 0.0;
                    }
                    return out;
                }

                const StepPropertyContainer resistivity;
            };
        } // ResFlowProfileControl
    } // Wells
} // GPN