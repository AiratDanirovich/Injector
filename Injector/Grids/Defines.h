#pragma once

#include <cassert>
#include <vector>

#include <Eigen/Core>

namespace GPN
{
    /// @brief It is convinient in debug mode to have a standard container.
    /// custom_vector type is introduced to overload operator() for element access via [].
    /// While release mode should be compiled with Eigen__Array as container.
    /// @tparam T value_type for std::vector<T>
    template<typename T>
    struct custom_vector : public std::vector<T>
    {        
        custom_vector(size_t size)
            : std::vector<T>(size)
        {}

        T& operator()(ptrdiff_t idx)
        { return (*this)[static_cast<size_t>(idx)];}
        
        const T& operator()(ptrdiff_t idx) const
        { return (*this)[static_cast<size_t>(idx)];}

        std::ptrdiff_t size() const
        {
            return static_cast<std::ptrdiff_t>(std::vector<T>::size());
        }
    };

    using RealType = double;
    using NodesContainer = custom_vector<RealType>; // use Eigen::ArrayX<RealType>; in Release
    using DualNodesContainer = NodesContainer;
    using MeshStepsContainer = NodesContainer;
    using DualStepsContainer = NodesContainer;
    using CellVolumeContainer = NodesContainer;
    using CellVolumeContainer2D = Eigen::ArrayXX<RealType>;


    struct Box
    {
        struct Left {RealType value;};
        struct Right {RealType value;};
        struct Top {RealType value;};
        struct Bottom {RealType value;};

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
        struct Step_R {RealType step;};
        struct Step_Z {RealType step;};

        Steps(Step_R step_r, Step_Z step_z)
            : step_r{step_r}, step_z{step_z}
        {
        }

        Step_R step_r;
        Step_Z step_z;
    };

    // x-axis goes up-down, South-North
    // y-axis goes left-right, East-West
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
            BCSouth(RealType fixed_x, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_x{fixed_x}
            {
            }

            RealType fixed_x;
        };

        struct BCNorth : public BoundaryCondition
        {
            BCNorth(RealType fixed_x, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_x{fixed_x}
            {
            }

            RealType fixed_x;
        };

        struct BCEast : public BoundaryCondition
        {
            BCEast(RealType fixed_y, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_y{fixed_y}
            {
            }

            RealType fixed_y;
        };

        struct BCWest : public BoundaryCondition
        {
            BCWest(RealType fixed_y, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_y{fixed_y}
            {
            }

            RealType fixed_y;
        };

        struct BCFunctorBase
        {
            virtual RealType operator()(RealType x, RealType y, RealType t) const = 0;
        };
    } // BoundaryConditions
} // EqSolver