#pragma once

#include <cassert>

#include <Eigen/Core>

namespace GPN
{
    using float_t = double;
    using NodesContainer = Eigen::ArrayX<float_t>;
    using MeshStepsContainer = NodesContainer;
    using CellVolumeContainer = NodesContainer;

    struct Box
    {
        struct Left {float_t value;};
        struct Right {float_t value;};
        struct Top {float_t value;};
        struct Bottom {float_t value;};

        Box(Left x_a, Right x_b, Top y_a, Bottom y_b)
            : x_a{x_a}, x_b{x_b}, y_a{y_a}, y_b{y_b}
        {
            assert(x_a.value < x_b.value);
            assert(y_a.value < y_b.value);
        };

        Left x_a; 
        Right x_b; 
        Top y_a; 
        Bottom y_b;
    };

    struct Steps
    {
        struct Step_R {float_t step;};
        struct Step_Z {float_t step;};

        Steps(Step_R step_r, Step_Z step_z)
            : step_r{step_r}, step_z{step_z}
        {
        }

        Step_R step_r;
        Step_Z step_z;
    };

    namespace BoundaryConditions
    {
        struct BoundaryCondition
        {
            enum BCType
            {
                first,
                second,
                third
            };
            BoundaryCondition(BCType type) : type{type} {}

            BCType type;
        };

        struct BCSouth : public BoundaryCondition
        {
            BCSouth(float_t fixed_x, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_x{fixed_x}
            {
            }

            float_t fixed_x;
        };

        struct BCNorth : public BoundaryCondition
        {
            BCNorth(float_t fixed_x, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_x{fixed_x}
            {
            }

            float_t fixed_x;
        };

        struct BCEast : public BoundaryCondition
        {
            BCEast(float_t fixed_y, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_y{fixed_y}
            {
            }

            float_t fixed_y;
        };

        struct BCWest : public BoundaryCondition
        {
            BCWest(float_t fixed_y, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_y{fixed_y}
            {
            }

            float_t fixed_y;
        };

        struct BCFunctorBase
        {
            virtual float_t operator()(float_t x, float_t y, float_t t) const = 0;
        };
    } // BoundaryConditions
} // EqSolver