#pragma once

#include <Injector/Grids/Defines.h>

#include <Injector/Solver/BoundaryConditions.hpp>

namespace GPN
{
    namespace Wells
    {
        namespace ResFlowProfileControl
        {
#pragma region BOUNDARY-CONDITION
            template <typename History_t, typename Grid2D_t>
            struct FunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
            {
                //    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
                FunctorBC(
                    const cptr<History_t> history,
                    const Logs::ExternalPressure &ext_pressure,
                    const StepPropertyContainer &rfp,
                    const cptr<const Grid2D_t> grid_ptr)
                    : history{history},
                      ext_pressure{ext_pressure},
                      grid_ptr{grid_ptr},
                      rfp{rfp}
                {
                    assert(rfp.sum() == 1.0);
                }

                RealType operator()(const ptrdiff_t z_id, const RealType r, const RealType,
                                    const BCType bc_type) const override
                {
                    if (r == grid_ptr->second_coord().dual_front())
                    {
                        assert(bc_type == BCType::second);
                        const auto out{rfp(z_id) * history->rate()};
                        return out;
                    }

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
                    // assume zero diffusion flux in hydrodynamic equation
                    assert(bc_type == BCType::second);
                    return 0.0;
                }

            protected:
                const Logs::ExternalPressure &ext_pressure;
                const cptr<const Grid2D_t> grid_ptr;
                const StepPropertyContainer rfp;
                const cptr<History_t> history;
            };

            struct RFPControlBC : public BoundaryConditions::GeneralBC
            {
                template <typename Grid2D_t>
                RFPControlBC(const cptr<Grid2D_t> &grid,
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
                typename Grid2D_t,
                typename Well_t>
            struct WellReservoirFlowProfileControl : public Well_t
            {
                using Well_t::RFP_weights;
                using Well_t::weights_sum;

                template <typename History_t>
                using functor_type = FunctorBC<History_t, Grid2D_t>;

                using hydro_bc_type = RFPControlBC;

                WellReservoirFlowProfileControl(const Well_t &well_base,
                                                const Properties::Rocks::RocksProps<Grid2D_t> &
                                                    rock_field_props,
                                                const cptr<Grid2D_t> grid2D_rocks)
                    : Well_t{well_base},
                      inv_mobility{set_inv_mobility(rock_field_props, well_base, grid2D_rocks)}
                {
                }

                StepPropertyContainer get_pressure_at_symmetry_axis(
                    const RealType rate,
                    const auto &ref_pressure) const
                {
                    return ref_pressure.col(0ll) + rate * inv_mobility;
                }

            private:
                static auto set_inv_mobility(
                    const Properties::Rocks::RocksProps<Grid2D_t> &rock_field_props,
                    const Well_t &well,
                    cptr<Grid2D_t> grid2D_rocks)
                {
                    const auto r3{grid2D_rocks->second_coord().mesh_nodes(0ll)};
                    const auto r2_face{grid2D_rocks->second_coord().dual_nodes(0ll)};
                    const auto two_pi{2.0 * std::numbers::pi_v<RealType>};
                    const auto is_permeable{rock_field_props.base_hydrodynamics.is_permeable.log_vals};

                    Eigen::ArrayX<RealType> out{well.RFP_weights / (two_pi / std::log(r3 / r2_face) *
                                                                    rock_field_props.mobility_axes2.col(Grid2D_t::l_margin) *
                                                                    grid2D_rocks->first_coord().volumes())};

                    for (auto row{0ll}; row < grid2D_rocks->first_coord().mesh_size(); ++row)
                    {
                        if (is_permeable(row) == 0.0)
                            out(row) = 0.0;
                    }
                    return out;
                }

                const StepPropertyContainer inv_mobility;
            };

        } // ResFlowProfileControl
    } // Wells
} // GPN