#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Solver/BoundaryConditions.hpp>

namespace GPN
{
    struct BCFunctor : public GPN::BoundaryConditions::GeneralBC::BCFunctorBase
    {
        BCFunctor(const RealType val) : out{BC_descriptor::BC_I(val)} {}

        BC_descriptor operator()(const ptrdiff_t, const RealType, const RealType,
                            const BCType) const override
        {
            return out;
        }

        BC_descriptor operator()(const RealType, const ptrdiff_t, const RealType,
                            const BCType) const override
        {
            return out;
        }

    protected:
        const BC_descriptor out;
    };
} // GPN