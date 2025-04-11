#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

#include <Eigen/Core>
// #include <Eigen/Dense>

#include "Defines.h"
#include "CoordinateTypes.h"

namespace GPN
{
    namespace Grids
    {
        /// @brief Struct defines the grid stencils, i.e.
        /// nodes that must be included as grid nodes.
        /// Further refinement is possible
        struct GridStencils
        {
            GridStencils(const std::vector<float_t> &nodes) noexcept 
                : mesh_nodes(nodes.size())
            {
                // at least two nodes required
                assert(nodes.size() > 1ull);
                // nodes must be strictly monotonuos
                for (size_t idx{0}; idx < nodes.size() - 1; ++idx)
                    assert(nodes[idx] < nodes[idx + 1]);

                // copy stencils to local container
                for (size_t idx{0}; idx < nodes.size(); ++idx)
                    mesh_nodes(idx) = nodes[idx];
            }

            GridStencils(GridStencils &&) noexcept = default;
            GridStencils(const GridStencils &) noexcept = default;
            GridStencils() = delete;

            const auto& get_nodes() const{
                return mesh_nodes; 
            }

        protected:
            // read-only
            NodesContainer mesh_nodes;
        };

        /// @brief Container for grid nodes.
        /// The nodes are somehow distributed, uniformly or non-uniformly.
        struct GridNodes
        {
            // Simple ctor without mesh refinement.
            // Only stencil nodes are used,
            // without mesh refinement
            GridNodes(const GridStencils &nodes) noexcept 
                : GridNodes{nodes, nodes.get_nodes()}
            {
            }

            template <typename RefinementPolicy>
            GridNodes(const GridStencils &nodes, RefinementPolicy &&refiner) noexcept
                : GridNodes{nodes, refiner.refine(nodes)}
            {
                // define refinement policy
                assert(false);
            }

            auto size() const { return mesh_nodes.size(); }

            GridNodes(GridNodes &&) noexcept = default;
            GridNodes(const GridNodes &) noexcept = default;
            GridNodes() = delete;

            const auto& get_nodes() const{
                return mesh_nodes; 
            }

            const auto operator[](size_t id) const{
                return mesh_nodes[id];
            }

        public:
            // mesh to be used in simulation
            NodesContainer mesh_nodes;
            // stencils of the mesh.
            // Here, jumps of physical properties occur
            GridStencils stencil_nodes;

        protected:
            GridNodes(const GridStencils &nodes, const NodesContainer& refined_mesh) noexcept 
            : stencil_nodes{nodes}
            , mesh_nodes{refined_mesh}
            {
                auto size{mesh_nodes.size()};
                // nodes must be strictly monotonuos
                for (auto idx{size-size}; idx < mesh_nodes.size() - 1; ++idx)
                    assert(mesh_nodes[idx] < mesh_nodes[idx + 1]);
            }
        };


        /// @brief Aggregates grid nodes and
        /// its associated properties, i.e.,
        /// volume per node, heat resistivity etc.
        template<typename CoordinateType_t>
        struct Grid1D
        {
            Grid1D(const GridNodes &nodes) noexcept 
            : mesh_nodes{nodes}
            , mesh_steps{CoordinateType_t::steps(nodes.mesh_nodes)}
            , cell_volumes{CoordinateType_t::volumes(nodes.mesh_nodes)}
            {
                assert(nodes.size() > 1ull);
            }

            Grid1D(Grid1D &&) noexcept = default;
            Grid1D(const Grid1D &) noexcept = default;
            Grid1D() = delete;

            auto node(auto id) const
            {
                assert(id < mesh_nodes.size());
                return mesh_nodes(id);
            }
            // step between nodes id and id+1
            auto step(auto id) const
            {
                assert(id < mesh_steps.size());
                return mesh_steps(id);
            }

            auto volume(auto idx) const{
                assert(idx < cell_volumes.size());
                return cell_volumes(idx);
            }

        protected:
            GridNodes mesh_nodes;
            MeshStepsContainer mesh_steps;
            CellVolumeContainer cell_volumes;
        };

    } // Grids
} // GPN