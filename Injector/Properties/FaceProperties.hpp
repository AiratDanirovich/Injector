#pragma once

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>

#include <Injector/Properties/PhysicalField.hpp>

namespace GPN
{
    namespace Logs
    {
        struct FaceInterpolator
        {
            template <typename Grid_t>
            static auto interpolate1D_z(
                const MeshNodesContainer &prop, const Grid_t &grid)
            {
                using CoordinateType_t = typename Grid_t::Axes;
                //    const auto &grid{log.grid};
                // if dual_size() == 2 --- no internal faces, only single cell exists
                assert(grid.dual_size() > 2ll);

                InternalFaceValues out(
                    grid.dual_size() - 2ll);

                for (auto id{0ll}; id < out.size(); ++id)
                    out(id) = CoordinateType_t::face_interpolator(
                        grid.mesh_nodes(id), grid.mesh_nodes(id + 1ll), grid.dual_nodes(id + 1ll), prop(id), prop(id + 1ll));

                return out;
            }

            template <typename Grid_t>
            static auto interpolate1D_r(
                const MeshNodesContainerT &prop, const Grid_t &grid)
            {
                using CoordinateType_t = typename Grid_t::Axes;
                // if dual_size() == 2 --- no internal faces, only single cell exists
                assert(grid.dual_size() > 2ll);

                InternalFaceValuesT out(
                    grid.dual_size() - 2ll);

                for (auto id{0ll}; id < out.size(); ++id)
                    out(id) = CoordinateType_t::face_interpolator(
                        grid.mesh_nodes(id), grid.mesh_nodes(id + 1ll), grid.dual_nodes(id + 1ll), prop(id), prop(id + 1ll));

                return out;
            }
        };

        // template <>
        // struct FaceInterpolator<CoordinateTypes::R_CylCoord>
        // {
        //     using Axes = CoordinateTypes::R_CylCoord;

        // public:
        //     template <typename Grid_t>
        //     static auto interpolate(
        //         const StepPropertyGrid &log,
        //         const Grid_t &grid)
        //     {
        //         // if dual_size() == 2 --- no internal faces
        //         assert(grid.dual_size() >= 2ll);

        //         FaceValuesContainer out(log.grid.mesh_size(), grid.dual_size() - 2ll);

        //         for (auto id{0ll}; id < out.cols(); ++id)
        //             out.col(id) = Axes::const_face_interpolator(
        //                 grid.mesh_nodes(id), grid.mesh_nodes(id + 1ll),
        //                 static_cast<StepPropertyContainer>(log.log_vals));
        //         return out;
        //     }
        // };

        // template <typename CoordinateType_t /* = CoordinateTypes::Z*/>
        // struct FaceInterpolatedProperty
        //     : public StepPropertyGrid,
        //       public FaceInterpolator<CoordinateType_t>
        // {
        //     FaceInterpolatedProperty(
        //         const StepPropertyGrid &vals) noexcept
        //         : StepPropertyGrid{vals},
        //           FaceInterpolator<CoordinateType_t>{vals}
        //     {
        //     }
        // };

        // using ZInterpolator =
        //     FaceInterpolatedProperty<CoordinateTypes::Z>;

        // using RInterpolator =
        //     FaceInterpolatedProperty<CoordinateTypes::R_CylCoord>;

    } // Logs

    namespace FaceProperties
    {
        struct FaceInterpolator
        {
            template <typename Grid2D_t>
            static auto interpolate2D_z(
                const CellNodesContainer2D &prop, const Grid2D_t &grid)
            {
                // if dual_size() == 2 --- no internal faces, only single cell exists
                assert(grid.coordinate<typename Grid2D_t::Axes1>().dual_size() > 2ll);

                FluxComponentContainer out(
                    grid.coordinate<typename Grid2D_t::Axes1>().dual_size() - 2ll,
                    grid.coordinate<typename Grid2D_t::Axes2>().mesh_size());

                for (auto id{0ll}; id < out.cols(); ++id)
                    out.col(id) = Logs::FaceInterpolator::interpolate1D_z(
                        prop.col(id),
                        grid.coordinate<typename Grid2D_t::Axes1>());

                return out;
            }

            template <typename Grid2D_t>
            static auto interpolate2D_r(
                const CellNodesContainer2D &prop, const Grid2D_t &grid)
            {
                // if dual_size() == 2 --- no internal faces, only single cell exists
                assert(grid.coordinate<typename Grid2D_t::Axes2>().dual_size() > 2ll);

                FluxComponentContainer out(
                    grid.coordinate<typename Grid2D_t::Axes1>().mesh_size(),
                    grid.coordinate<typename Grid2D_t::Axes2>().dual_size() - 2ll);

                for (auto id{0ll}; id < out.rows(); ++id)
                // take a row of mesh-node values
                // and interpolate at faces to get dual-node values
                    out.row(id) = Logs::FaceInterpolator::interpolate1D_r(
                        prop.row(id),
                        grid.coordinate<typename Grid2D_t::Axes2>());

                return out;
            }
        };

        template <typename Grid_t>
        struct FaceInterpolatedField
        {
            FaceInterpolatedField(
                const FaceValuesContainer &face_vals_axes1,
                const FaceValuesContainer &face_vals_axes2)
                : face_vals_axes1{face_vals_axes1},
                  face_vals_axes2{face_vals_axes2}
            {
                assert(face_vals_axes2.rows() > 0ll);
                assert(face_vals_axes1.cols() > 0ll);
            }

            // interpolated values at faces normal to Axes1
            FaceValuesContainer face_vals_axes1;
            // interpolated values at faces normal to Axes2
            FaceValuesContainer face_vals_axes2;
        };

        struct FaceInterpolatedFieldFactory
        {
            template <typename Grid_t>
            static auto create(
                const Properties::Field<Grid_t> &property,
                const cptr<Grid_t> grid)
            {
                // FaceValuesContainer face_vals_axes1{
                //     FaceInterpolator::interpolate2D_z(
                //         property.values(),
                //         *grid)};

                // FaceValuesContainer face_vals_axes2{
                //     FaceInterpolator::interpolate2D_r(
                //         property.values(),
                //         *grid)};

                return FaceInterpolatedField<Grid_t>{
                    FaceValuesContainer{
                        FaceInterpolator::interpolate2D_z(
                            property.values(),
                            *grid)},
                    FaceValuesContainer{
                        FaceInterpolator::interpolate2D_r(
                            property.values(),
                            *grid)}};
            }
        };

        template <typename Grid_t>
        using HeatConductivity = FaceInterpolatedField<Grid_t>;
        
        template <typename Grid_t>
        using MediumHeatConductivity = HeatConductivity<Grid_t>;

        // struct ThermalDiffusivity
        //     : public Properties::Field<
        //           Grids::StructuredCylinderGrid2DAxisymmetric>
        // {
        //     ThermalDiffusivity(const Properties::MediumHeatVolumetricCapacity &capacity,
        //                        const Properties::HeatConductivity &conductivity)
        //         : Field{
        //               conductivity / capacity}
        //     {
        //     }
        // };

    } // FaceProperties

} // GPN