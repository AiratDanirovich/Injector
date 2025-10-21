#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/History/History.hpp>

#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/State2D.hpp>

#include <Injector/Properties/Logs.hpp>

using VR = std::vector<GPN::RealType>;

namespace GPN
{
  /// @brief Initial temperature is assumed to be constant
  struct FunctorIC : public GPN::InitialConditions::ICFunctorBase
  {
    FunctorIC(const Logs::Geotherma &geotherma)
        : geotherma{geotherma}
    {
    }
    RealType operator()(const ptrdiff_t z_id, const ptrdiff_t, GPN::RealType) const override
    {
      return geotherma(z_id);
    }

  protected:
    const Logs::Geotherma &geotherma;
  };

  template <typename Grid_t_ptr>
  auto ICFactory(RealType t0, const Grid_t_ptr grid, const Logs::Geotherma &geotherma)
  {
    return EqSolver::State::State2D{EqSolver::State::State2D::FillWithFunctor(*grid, FunctorIC{geotherma}, t0)};
  }

  template <typename Well_t, typename Hydro_t>
  struct FunctorBC : public BoundaryConditions::GeneralBC::BCFunctorBase
  {
    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
    using BCType = BoundaryConditions::GeneralBC::BoundaryCondition::BCType;
    using ConvectionFieldFactory_t =
        FaceProperties::IncompressibleRatesFactory<
            Grid2D_t, Well_t, History, PhasePropertiesJT, Hydro_t>;
    FunctorBC(
        const ptr<ConvectionFieldFactory_t> flow_field, // volumetric heat flow rate
        const Logs::Geotherma &geotherma,
        const cptr<const Grid2D_t> grid_ptr)
        : flow_field{flow_field},
          geotherma{geotherma},
          grid_ptr{grid_ptr}
    {
    }

    RealType operator()(const ptrdiff_t z_id, const RealType r, const RealType t,
                        const BCType bc_type) const override
    {
      if (r == grid_ptr->second_coord().dual_front())
        return 0.0; // bc at the axis of symmetry, r == 0.0

      if (r == grid_ptr->second_coord().dual_back())
      {                // bc at the external contour
        if (bc_type == // geotherma is set for producer
            BCType::first)
          return geotherma(z_id);
        else if (bc_type == // zero diffusion flux for injector
                 BCType::second)
          return 0.0;
      }

      assert(false);
      return 0.0;
    }

    RealType operator()(const RealType z, const ptrdiff_t r_id, const RealType t,
                        const BCType) const override
    {
      if (z == grid_ptr->first_coord().dual_front())
      { // inflow with temperature from history,
        // outflow is accounted for in the matrix
        return std::max(0.0, flow_field->get_heat_flow_in_axes1()(0ll, r_id)) * flow_field->get_temperature();
      }

      if (z == grid_ptr->first_coord().dual_back())
      {
        // outflow -- duffusion flux is zero, min -> 0.0
        // inflow -- geotherm inflows from the bottom hole
        return std::min(0.0, flow_field->get_heat_flow_in_axes1()(
                                 grid_ptr->first_coord().dual_size() - 1ll, r_id)) *
               geotherma.log_vals.tail(1ll)(0ll);
      }

      assert(false);
      return 0.0;
    }

  protected:
    const Logs::Geotherma &geotherma;
    const cptr<const Grid2D_t> grid_ptr;
    const ptr<ConvectionFieldFactory_t> flow_field;
  };

} // GPN