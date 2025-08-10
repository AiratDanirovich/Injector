#pragma once

#include <vector>
#include <cassert>

#include <Injector/Grids/Defines.h>
#include <Injector/Model/Completion.hpp>

namespace GPN
{
    namespace Completion
    {
        struct ExtrudedRingFactory
        {
            template <typename Grid_t>
            static ExtrudedRing create_ring(
                const FlowRing<Ring> &flo_ring,
                const TubeRing<Ring> &tub_ring,
                const ColumnRing<Ring> &col_ring,
                const Grid_t &grid_z)
            {
                using namespace Eigen;

                const auto size{grid_z.mesh_size()};
                ArrayX<RealType> thickness(size), inner_radius(size);
                inner_radius = 0.0;
                for (auto id{0ll}; id < size; ++id)
                {
                    thickness(id) = grid_z.mesh_nodes(id) < tub_ring.max_depth
                                        ? flo_ring.thickness
                                        : col_ring.inner_radius;
                }

                return ExtrudedRing{flo_ring, inner_radius, thickness};
            }

            template <typename Grid_t, typename Ring_t>
            static ExtrudedRing create_ring(
                const TubeRing<Ring_t> &tub_ring,
                const ColumnRing<Ring_t> &col_ring,
                const Grid_t &grid_z)
            {
                using namespace Eigen;

                const auto size{grid_z.mesh_size()};
                ArrayX<RealType> thickness(size), inner_radius(size);
                for (auto id{0ll}; id < size; ++id)
                {
                    thickness(id) = grid_z.mesh_nodes(id) < tub_ring.depth
                                        ? tub_ring.thickness
                                        : 0.0;
                    inner_radius(id) = grid_z.mesh_nodes(id) < tub_ring.depth
                                           ? tub_ring.inner_radius
                                           : col_ring.inner_radius;
                }

                return ExtrudedRing{tub_ring, inner_radius, thickness};
            }

            template <typename Grid_t, typename Ring_t>
            static ExtrudedRing create_ring(
                const AnnulusRing<Ring_t> &ann_ring,
                const TubeRing<Ring_t> &tub_ring,
                const ColumnRing<Ring_t> &col_ring,
                const Grid_t &grid_z)
            {
                using namespace Eigen;

                const auto size{grid_z.mesh_size()};
                ArrayX<RealType> thickness(size), inner_radius(size);
                for (auto id{0ll}; id < size; ++id)
                {
                    thickness(id) = grid_z.mesh_nodes(id) < tub_ring.depth
                                        ? ann_ring.thickness
                                        : 0.0;
                    inner_radius(id) = grid_z.mesh_nodes(id) < tub_ring.depth
                                           ? ann_ring.inner_radius
                                           : col_ring.inner_radius;
                }

                return ExtrudedRing{ann_ring, inner_radius, thickness};
            }

            template <typename Grid_t, typename Ring_t>
            static ExtrudedRing create_ring(
                const ColumnRing<Ring_t> &ring,
                const Grid_t &grid_z)
            {
                using namespace Eigen;

                const auto size{grid_z.mesh_size()};
                const auto thickness{ArrayX<RealType>::Constant(size, ring.thickness)};
                const auto inner_radius{ArrayX<RealType>::Constant(size, ring.inner_radius)};

                return ExtrudedRing{ring, inner_radius, thickness};
            }

            template <typename Grid_t, typename Ring_t>
            static ExtrudedRing create_ring(
                const CementRing<Ring_t> &ring,
                const Grid_t &grid_z)
            {
                using namespace Eigen;

                const auto size{grid_z.mesh_size()};
                const auto thickness{ArrayX<RealType>::Constant(size, ring.thickness)};
                const auto inner_radius{ArrayX<RealType>::Constant(size, ring.inner_radius)};

                return ExtrudedRing{ring, inner_radius, thickness};
            }
        };

        struct ExtrudedCasingFactory
        {
            template <typename Grid_t>
            static ExtrudedCasing create(
                const Casing<Ring> &casing,
                const Grid_t &grid_z)
            {
                std::vector<ExtrudedRing> extruded_casing;
                extruded_casing.reserve(MaterialType::Size);

                extruded_casing.emplace_back(
                    ExtrudedRingFactory::create_ring(
                        FlowRing{casing[MaterialType::Flow]},
                        TubeRing{casing[MaterialType::Tube]},
                        ColumnRing{casing[MaterialType::Column]},
                        grid_z));

                extruded_casing.emplace_back(
                    ExtrudedRingFactory::create_ring(
                        TubeRing{casing[MaterialType::Tube]},
                        ColumnRing{casing[MaterialType::Column]},
                        grid_z));

                extruded_casing.emplace_back(
                    ExtrudedRingFactory::create_ring(
                        AnnulusRing{casing[MaterialType::Annulus]},
                        TubeRing{casing[MaterialType::Tube]},
                        ColumnRing{casing[MaterialType::Column]},
                        grid_z));

                extruded_casing.emplace_back(
                    ExtrudedRingFactory::create_ring(
                        ColumnRing{casing[MaterialType::Column]},
                        grid_z));

                extruded_casing.emplace_back(
                    ExtrudedRingFactory::create_ring(
                        CementRing{casing[MaterialType::Cement]},
                        grid_z));

                for (auto id{0ll}; id < grid_z.mesh_size(); ++id)
                {
                    for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Cement; ++i)
                    {
                        assert(is_tight_casing(extruded_casing[i].inner_radius(id), extruded_casing[i - 1ll].outer_radius(id)));
                    }
                    for (const auto &r : extruded_casing)
                    {
                        assert(std::abs(r.inner_radius(id) + r.thickness(id) - r.outer_radius(id)) < 1e-12);
                    }
                }

                return ExtrudedCasing{extruded_casing};
            }

            // template <typename Grid_t>
            // static ExtrudedCasing create(
            //     const Casing<VarRing> &casing,
            //     const Grid_t &grid_z)
            // {
            //     std::vector<ExtrudedRing> extruded_casing;
            //     extruded_casing.reserve(MaterialType::Size);

            //     extruded_casing.emplace_back(
            //         ExtrudedRingFactory::create_ring(
            //             FlowRing{casing[MaterialType::Flow]},
            //             TubeRing{casing[MaterialType::Tube]},
            //             ColumnRing{casing[MaterialType::Column]},
            //             grid_z));

            //     extruded_casing.emplace_back(
            //         ExtrudedRingFactory::create_ring(
            //             TubeRing{casing[MaterialType::Tube]},
            //             ColumnRing{casing[MaterialType::Column]},
            //             grid_z));

            //     extruded_casing.emplace_back(
            //         ExtrudedRingFactory::create_ring(
            //             AnnulusRing{casing[MaterialType::Annulus]},
            //             TubeRing{casing[MaterialType::Tube]},
            //             ColumnRing{casing[MaterialType::Column]},
            //             grid_z));

            //     extruded_casing.emplace_back(
            //         ExtrudedRingFactory::create_ring(
            //             ColumnRing{casing[MaterialType::Column]},
            //             grid_z));

            //     extruded_casing.emplace_back(
            //         ExtrudedRingFactory::create_ring(
            //             CementRing{casing[MaterialType::Cement]},
            //             grid_z));

            //     for (auto id{0ll}; id < grid_z.mesh_size(); ++id)
            //     {
            //         for (ptrdiff_t i{MaterialType::Tube}; i <= MaterialType::Cement; ++i)
            //         {
            //             assert(is_tight_casing(extruded_casing[i].inner_radius(id), extruded_casing[i - 1ll].outer_radius(id)));
            //         }
            //         for (const auto &r : extruded_casing)
            //         {
            //             assert(std::abs(r.inner_radius(id) + r.thickness(id) - r.outer_radius(id)) < 1e-12);
            //         }
            //     }

            //     return ExtrudedCasing{extruded_casing};
            // }
        };

    } // Completion
} // GPN