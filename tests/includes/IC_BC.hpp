#pragma once

#include <Injector/Grids/Defines.h>

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
  struct FunctorBC : public BoundaryConditions::BCFunctorBase
  {
    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
    using ConvectionFieldFactory_t =
        FaceProperties::IncompressibleRatesFactory<
            Grid2D_t, Well_t, PhasePropertiesJT, Hydro_t>;
    FunctorBC(
        const ConvectionFieldFactory_t &flow_field, // volumetric heat flow rate
        const cptr<const Grid2D_t> grid_ptr)
        : flow_field{flow_field},
          grid_ptr{grid_ptr}
    {
    }

    RealType operator()(const ptrdiff_t z_id, const RealType r, const RealType t,
                        const BoundaryConditions::BoundaryCondition::BCType bc_type =
                            BoundaryConditions::BoundaryCondition::BCType::second) const override
    {
      if (r == grid_ptr->second_coord().dual_front())
        return flow_field.get_heat_flow_in_axes2()(z_id, 0ll) * flow_field.get_temperature();

      if (r == grid_ptr->second_coord().dual_back())
      {
        if (bc_type ==
            BoundaryConditions::BoundaryCondition::BCType::first)
          return 0.0;
        else if (bc_type ==
                 BoundaryConditions::BoundaryCondition::BCType::second)
          return 0.0;
        else
          return 0.0;
      }

      assert(false);
      return 0.0;
    }

    RealType operator()(const RealType z, const ptrdiff_t r_id, const RealType t,
                        const BoundaryConditions::BoundaryCondition::BCType bc_type =
                            BoundaryConditions::BoundaryCondition::BCType::second) const override
    {
      if (z == grid_ptr->first_coord().dual_front())
        return flow_field.get_heat_flow_in_axes1()(0ll, r_id) * flow_field.get_temperature();

      if (z == grid_ptr->first_coord().dual_back())
        return 0.0;

      assert(false);
      return 0.0;
    }

  protected:
    const cptr<const Grid2D_t> grid_ptr;
    const ConvectionFieldFactory_t &flow_field;
  };

} // GPN