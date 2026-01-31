#pragma once

#include <cassert>

#include <Injector/Grids/Defines.h>

#include <Injector/Model/Well/DefaultWellNmerics.hpp>

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
            template <typename History_t, typename Grid2D_t>
            struct BotHolePresFunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
            {
                /// @brief 
                /// @param history 
                /// @param ext_pressure 
                /// @param PI productivity index that multiplies pressure difference at sandface
                /// @param grid_ptr 
                BotHolePresFunctorBC(
                    const ptr<const History_t> history,
                    const Logs::ExternalPressure &ext_pressure,
                    const StepPropertyContainer& PI,
                    const ptr<const Grid2D_t> grid_ptr)
                    : BoundaryConditions::GeneralBC::BCFunctorBase{},
                      history{history},
                      ext_pressure{ext_pressure},
                      PI{PI},
                      grid_ptr{grid_ptr}
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
                        const auto out{history->pressure()};
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

            protected:
                const Logs::ExternalPressure &ext_pressure;
                const StepPropertyContainer& PI;
                const ptr<const Grid2D_t> grid_ptr;
                const ptr<const History_t> history;
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
                typename Grid2D_t,
                typename CrossFlow_t>
            struct WellBottomHolePressureControl : 
                public CrossFlow_t, 
                public DefaultWellNumerics<0ll>
            {
                using functor_type = BotHolePresFunctorBC<History_t, Grid2D_t>;
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
                    const cptr<Grid2D_t> grid2D_rocks)
                    : CrossFlow_t{well_base},
                      size{grid2D_rocks->first_coord().mesh_size()},
                      PI{set_productivity_index(rock_field_props, grid2D_rocks)},
                      flow_axes1_value{
                          FaceValuesContainer::Zero(
                              grid2D_rocks->first_coord().mesh_size() + 1ll,
                              Grid2D_t::l_margin)},
                      flow_axes2_value{
                          FaceValuesContainer::Zero(
                              grid2D_rocks->first_coord().mesh_size(),
                              Grid2D_t::l_margin + 1ll)},
                      history{history},
                      rock_field_props{rock_field_props},
                      grid2D_rocks{grid2D_rocks},
                      is_permeable{rock_field_props.base_hydrodynamics.is_permeable}
                {
                    const_cast<ptr<const hydro_bc_type>&>(hydro_bc) = 
                          std::make_shared<const hydro_bc_type>(
                              grid2D_rocks,
                              std::make_shared<const functor_type>(
                                  history, rock_field_props.base_hydrodynamics.ext_pressure,
                                  PI,
                                  grid2D_rocks));
                }

                template <typename HistoryRecord_t>
                const auto get_RFP(const HistoryRecord_t &record) const
                {
                    return this->rfp;
                }
                template <typename HistoryRecord_t>
                const auto get_WFP(const HistoryRecord_t &record) const
                {
                    return this->wfp;
                }
                template <typename HistoryRecord_t>
                const auto get_verticle_cement_flow(const HistoryRecord_t &record) const
                {
                    return this->verticle_flux_in_cement;
                }
                template <typename HistoryRecord_t>
                const auto get_verticle_well_flow(const HistoryRecord_t &record) const
                {
                    return this->verticle_flux_in_well;
                }

                template <typename HistoryRecord_t>
                StepPropertyContainer get_pressure_at_sandface(
                    const HistoryRecord_t &record,
                    const auto &) const
                {
                    // Pressure at the level of NON-permeable layers is assumed zero.
                    // Pressure at permeable layers is equal to history->pressure()
                    return StepPropertyContainer::Constant(size, record.pressure) * is_permeable.log_vals;
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
                    this->set_flux(set_RFP(record, collector_pressure));
                    // set verticle flux
                    flow_axes1_value.col(0ll) = get_verticle_well_flow(record);
                    flow_axes1_value.col(1ll) = 0.0;
                    flow_axes1_value.col(2ll) = get_verticle_cement_flow(record);
                    // set radial flux
                    static_assert(Grid2D_t::l_margin == 3ll);

                    flow_axes2_value.col(0ll) = 0.0;
                    flow_axes2_value.col(1ll) = get_WFP(record);
                    flow_axes2_value.col(2ll) = get_WFP(record);
                    flow_axes2_value.col(3ll) = get_RFP(record);
                }

                FaceValuesContainer flow_axes1_value, flow_axes2_value;
                const cptr<Grid2D_t> grid2D_rocks;
                const cptr<History_t> history;
                const Properties::Rocks::RocksProps<Grid2D_t> &
                    rock_field_props;

            protected:
                template <typename HistoryRecord_t>
                const auto set_RFP(
                    const HistoryRecord_t &record,
                    const auto &collector_pressure) const
                {
                    const StepPropertyContainer depression_at_sandface{
                        -(collector_pressure.col(0ll) - record.pressure*is_permeable.log_vals)};
                    assert(
                        std::all_of(
                            depression_at_sandface.cbegin(), 
                            depression_at_sandface.cend(), 
                            [](const RealType v){
                                return !std::isnan(v);
                            }));

                    return Logs::RFPFactory::create_from_container<Logs::RFP>(
                        (PI * (collector_pressure.col(0ll) - record.pressure)).eval(),
                        is_permeable);
                }

                static auto set_productivity_index(
                    const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                    const ptr<const Grid2D_t> grid2D_rocks)
                {
                    // center of cell next to the well
                    const auto r3{grid2D_rocks->second_coord().mesh_nodes(0ll)};
                    // sandface
                    const auto r2_face{grid2D_rocks->second_coord().dual_nodes(0ll)};
                    const auto two_pi{2.0 * std::numbers::pi_v<RealType>};
                    const auto is_permeable{rock_field_props.base_hydrodynamics.is_permeable.log_vals};

                    Eigen::ArrayX<RealType> out{
                        two_pi / std::log(r3 / r2_face) *
                        rock_field_props.mobility_axes2.col(Grid2D_t::l_margin) *
                        grid2D_rocks->first_coord().volumes()};

                    for (auto row{0ll}; row < grid2D_rocks->first_coord().mesh_size(); ++row)
                    {
                        if (is_permeable(row) == 0.0)
                            out(row) = 0.0;
                    }
                    return out;
                }



                const std::ptrdiff_t size;
                const StepPropertyContainer PI;
                const Logs::IsPermeable &is_permeable;
            };
        } // BotHolePresControl
    } // Wells
} // GPN