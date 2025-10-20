#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Solver/BoundaryConditions.hpp>

namespace GPN
{
    struct BCFunctor : public GPN::BoundaryConditions::GeneralBC::BCFunctorBase
    {
        BCFunctor(const RealType val) : val{val} {}

        RealType operator()(const ptrdiff_t, const RealType, const RealType,
                            const BCType) const override
        {
            return val;
        }

        RealType operator()(const RealType, const ptrdiff_t, const RealType,
                            const BCType) const override
        {
            return val;
        }

    protected:
        const RealType val;
    };
} // GPN