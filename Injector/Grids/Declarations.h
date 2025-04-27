#pragma once

// https://blog.andreiavram.ro/object-has-method-cpp20-concepts/

#include <concepts>

using RealType = double;

namespace GPN
{
    template <typename Axes1_t, typename Axes2_t>
    struct CoordinateSystem2D;

    template <typename Axes1_t, typename Axes2_t, typename Axes3_t>
    struct CoordinateSystem3D;

    namespace CoordinateTypes
    {
        template <typename Axes_t>
        concept ICoordinate = requires(
            Axes_t, const DualNodesContainer &nodes,
            RealType xL, RealType xR, RealType xMid,
            RealType valL, RealType valR,
            const Eigen::ArrayX<RealType> &val) {
            { Axes_t::control_volumes(nodes) } -> std::convertible_to<ControlVolumesContainer>;
            { Axes_t::face_interpolator(
                xL, xR, xMid, valL, valR) } -> std::same_as<RealType>;
            { Axes_t::const_face_interpolator(
                xL, xR, val) } -> std::convertible_to<Eigen::ArrayX<RealType>>;
            { Axes_t::dual_steps(nodes) } -> std::convertible_to<DualStepsContainer>;
            { Axes_t::cell_centers(nodes) } -> std::convertible_to<MeshNodesContainer>;
            { Axes_t::mesh_steps(nodes) } -> std::convertible_to<MeshStepsContainer>;
        };

        struct GeneralCoordinate;

        struct CartesianCoordinate;
        struct X;
        struct Y;
        struct Z;
        struct R_CylCoord;

        struct Time;
    } // CoordinateTypes

    struct Cartesian2DCoordinates;

    struct CylinderCoordinates;

    struct Cartesian3DCoordinates;

    namespace Grids // 1D space
    {
        struct GridDualStencils;
        struct TemporalGridDualStencils;
        struct GridDual;
        struct TemporalGridDual;

        template <typename CoordinateType_t>
            requires CoordinateTypes::ICoordinate<CoordinateType_t>
        struct AxesGrid;

        struct RGrid;
        struct ZGrid;
        struct XGrid;
        struct YGrid;
    } // Grids

    namespace Grids // 2D space
    {
        template <typename CoordinateSystem_t>
        concept IStructuredGrid2D = requires {
            typename CoordinateSystem_t::Axes1;
            typename CoordinateSystem_t::Axes2;
            CoordinateSystem_t::Dim();
        };

        template <typename CoordinateSystem_t>
            requires IStructuredGrid2D<CoordinateSystem_t>
        struct StructuredGrid2D;

        struct StructuredCylinderGrid2DAxisymmetric;
        struct StructuredXYGrid2D;
    } // Grids
} // GPN