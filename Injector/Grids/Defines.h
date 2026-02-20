#pragma once

#include <cassert>
#include <vector>
#include <memory>

#include <Eigen/Core>

namespace GPN
{
    /// @brief It is convinient in debug mode to have a standard container.
    /// custom_vector type is introduced to overload operator() for element access via [].
    /// While release mode should be compiled with Eigen__Array as container.
    /// @tparam T value_type for std::vector<T>
    template <typename T>
    struct custom_vector : public std::vector<T>
    {
        custom_vector(size_t size)
            : std::vector<T>(size)
        {
        }

        T &operator()(ptrdiff_t idx)
        {
            return (*this)[static_cast<size_t>(idx)];
        }

        const T &operator()(ptrdiff_t idx) const
        {
            return (*this)[static_cast<size_t>(idx)];
        }

        std::ptrdiff_t size() const
        {
            return static_cast<std::ptrdiff_t>(std::vector<T>::size());
        }

        custom_vector operator*(T a) const
        {
            custom_vector out(size());
            for (std::ptrdiff_t id{0}; id < size(); ++id)
                out(id) = (*this)(id)*a;

            return out;
        }

        operator Eigen::ArrayX<T>() const
        {
            Eigen::ArrayX<T> out(this->size());
            for (ptrdiff_t id{0}; id < out.size(); ++id)
                out(id) = (*this)(id);
            return out;
        }
    };

    using RealType = double;
    using MeshNodesContainer = Eigen::ArrayX<RealType>;        // column array
    using MeshNodesContainerT = Eigen::Array<RealType, 1, -1>; // row array
    using LogValuesContainer = MeshNodesContainer;

    using StepPropertyContainer = Eigen::ArrayX<RealType>;
    using InternalFaceValues = Eigen::ArrayX<RealType>;
    using InternalFaceValuesT = Eigen::Array<RealType, 1, -1>; // row array

    /// @brief Container for the dual nodes coordinates
    struct DualNodesContainer : public MeshNodesContainer
    {
        using Base = MeshNodesContainer;
        using MeshNodesContainer::MeshNodesContainer;
    };
    using MeshStepsContainer = MeshNodesContainer;

    /// @brief Container for the normal distance
    /// between two faces of control volume
    struct DualStepsContainer : public MeshNodesContainer
    {
        using MeshNodesContainer::MeshNodesContainer;
        DualStepsContainer(const std::vector<RealType> &data)
            : MeshNodesContainer(data.size())
        {
            for (std::ptrdiff_t i{0}; i < static_cast<ptrdiff_t>(data.size()); ++i)
                (*this)[i] = data[i];
        }
    };
    using ControlVolumesContainer = MeshNodesContainer;

    using CellVolumeContainer2D = Eigen::ArrayXX<RealType>;
    using CellNodesContainer2D = CellVolumeContainer2D;
    using FluxComponentContainer = Eigen::ArrayXX<RealType>;

    using GridNodeValues2D = Eigen::ArrayXX<RealType>;
    using FaceValuesContainer = GridNodeValues2D;


    template <typename T>
    using ptr = std::shared_ptr<T>;

    template <typename T>
    using cptr = std::shared_ptr<T>;

    struct Gravity
    {
        constexpr static RealType value() noexcept
        {
            return (RealType)9.81;
        }
    };

    //    struct Directions
    //    {
    enum struct Directions
    {
        x1,
        x2,
        size
    };
    //    };

    struct ScalarParameter
    {
        RealType operator()(std::ptrdiff_t)
        {
            return value;
        }

        RealType operator()(RealType)
        {
            return value;
        }
        operator RealType() const { return value; }

    protected:
        RealType value;
    };

    struct LogParameter
    {
        LogParameter(const MeshNodesContainer &values) noexcept
            : values{values}
        {
        }

        RealType operator()(std::ptrdiff_t idx)
        {
            assert(idx < values.size());
            return values(idx);
        }

        // RealType operator()(RealType coord)
        // {
        //     return values(idx);
        // }

        operator MeshNodesContainer() const { return values; }

    protected:
        MeshNodesContainer values;
    };

    struct Segment
    {
        RealType start, end;
        Segment(RealType start, RealType end)
            : start{start},
              end{end}
        {
            assert(start < end);
        }
    };

    struct Box
    {
        Segment axes1, axes2;

        Box(const Segment &axes1, const Segment &axes2)
            : axes1{axes1}, axes2{axes2}
        {
        }
    };

    struct Steps
    {
        struct Step_R
        {
            RealType step;
        };
        struct Step_Z
        {
            RealType step;
        };

        Steps(Step_R step_r, Step_Z step_z)
            : step_r{step_r}, step_z{step_z}
        {
        }

        Step_R step_r;
        Step_Z step_z;
    };

    namespace InitialConditions
    {
        struct ICFunctorBase
        {
            virtual RealType operator()(const ptrdiff_t, const ptrdiff_t, const RealType) const = 0;
        };
        struct ICFunctorBase1D
        {
            virtual RealType operator()(const ptrdiff_t, const RealType) const = 0;
        };
    } // InitialConditions
} // EqSolver