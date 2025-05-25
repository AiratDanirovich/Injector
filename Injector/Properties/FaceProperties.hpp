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
                    out.col(id) = interpolate1D_z(
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
                    out.row(id) = interpolate1D_z(
                        prop.row(id),
                        grid.coordinate<typename Grid2D_t::Axes2>());

                return out;
            }

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

    //     namespace FaceProperties
    //     {
    //         template <typename Grid_t = Grids::StructuredCylinderGrid2DAxisymmetric>
    //         struct FaceInterpolatedField : public Properties::Field<Grid_t>
    //         {
    //             FaceInterpolatedField(
    //                 const Properties::Field<Grid_t> &property,
    //                 const cptr<Grid_t> grid)
    //                 : face_vals_axes1(
    //                       property.face_values.size(),
    //                       grid->second_coord.mesh_size()),
    //                   // #pragma region AXES2-FACEVALUES
    //                   face_vals_axes2{
    //                       Logs::FaceInterpolator<
    //                           typename Grid_t::Axes2>::interpolate(property,
    //                                                                grid->second_coord)}
    //             // #pragma endregion
    //             {
    // #pragma region AXES1-FACEVALUES
    //                 assert(property.face_values.size() >= 0ll);
    //                 // extrapolate Axes1-facevals as const value in the Axes2 direction
    //                 for (std::ptrdiff_t col{0ll}; col < face_vals_axes1.cols(); ++col)
    //                     face_vals_axes1.col(col) = property.face_values;

    //                 assert(face_vals_axes1.cols() > 0ll);
    //                 for (std::ptrdiff_t row{0ll}; row < face_vals_axes1.rows(); ++row)
    //                     for (std::ptrdiff_t col{1ll}; col < face_vals_axes1.cols(); ++col)
    //                         assert(face_vals_axes1(row, 0) == face_vals_axes1(row, col));
    // #pragma endregion
    //             }

    //             using typename Properties::Field<Grid_t>::Grid_type;

    //             // interpolated values at faces normal to Axes1
    //             FaceValuesContainer face_vals_axes1;
    //             // interpolated values at faces normal to Axes2
    //             FaceValuesContainer face_vals_axes2;
    //         };

    //         struct HeatConductivity
    //             : public FaceInterpolatedField<
    //                   Grids::StructuredCylinderGrid2DAxisymmetric>
    //         {
    //             using FaceInterpolatedField<
    //                 Grids::StructuredCylinderGrid2DAxisymmetric>::FaceInterpolatedField;
    //         };

    //         struct ThermalDiffusivity
    //             : public Properties::Field<
    //                   Grids::StructuredCylinderGrid2DAxisymmetric>
    //         {
    //             ThermalDiffusivity(const Properties::MediumHeatVolumetricCapacity &capacity,
    //                                const Properties::HeatConductivity &conductivity)
    //                 : Field{
    //                       conductivity / capacity}
    //             {
    //             }
    //         };

    //     } // FaceProperties

} // GPN