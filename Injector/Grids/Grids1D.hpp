#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

#include <Eigen/Core>
// #include <Eigen/Dense>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateTypes.h>

namespace GPN
{
    namespace Grids
    {
        /// @brief Container for the dual grid stencils, i.e.
        /// nodes that must be included as grid nodes in dual grid,
        /// i.e., the cell boundaries.
        /// Physically, they are the boundaries of various layers.
        /// Further refinement is possible.
        struct GridDualStencils
        {
            GridDualStencils(const std::vector<RealType> &nodes) noexcept
                : dual_nodes(nodes.size())
            {
#pragma region ASSERTIONS
                // at least two nodes are required,
                // min and max boundares of computation domain
                // in that particular direction
                assert(nodes.size() > 1ull);
                // stencils must be strictly monotonuos
                for (size_t idx{0}; idx < nodes.size() - 1; ++idx)
                    assert(nodes[idx] < nodes[idx + 1]);
#pragma endregion

                // copy dual mesh stencils to local container
                for (size_t idx{0ull}; idx < nodes.size(); ++idx)
                    dual_nodes(idx) = nodes[idx];
            }

            GridDualStencils(GridDualStencils &&) noexcept = default;
            GridDualStencils(const GridDualStencils &) noexcept = default;
            GridDualStencils() = delete;

            const auto &get_dual_nodes() const
            {
                return dual_nodes;
            }

            auto operator()(auto idx) const
            {
                assert(static_cast<ptrdiff_t>(idx) < static_cast<ptrdiff_t>(dual_nodes.size()));
                return dual_nodes(idx);
            }

            auto size() const
            {
                return dual_nodes.size();
            }

        protected:
            // read-only
            DualNodesContainer dual_nodes;
        };

        /// @brief Container for dual grid nodes (refined as well as stencils).
        /// The nodes are somehow distributed, 
        /// uniformly or non-uniformly, between adjuscent stencils.
        struct GridDual
        {
        private:
            using CoordinateType = CoordinateTypes::GeneralCoordinate;
        public:
            /// @brief Simple ctor without mesh refinement.
            /// Only stencil nodes are used,
            /// without mesh refinement
            /// @param nodes stencils of dual mesh
            GridDual(const GridDualStencils &nodes) noexcept
                : GridDual{
                      nodes,                  // stencils for the dual mesh
                      nodes.get_dual_nodes()} // nodes of the dual mesh --- same as stencils
            {
            }

            /// @brief Ctor with mesh refinement
            /// @tparam RefinementPolicy Type of mesh refinement policy
            /// @param nodes
            /// @param policy
            template <typename RefinementPolicy>
            GridDual(
                const GridDualStencils &nodes,
                RefinementPolicy &&policy) noexcept
                : GridDual{
                      nodes,               // nodes of dual mesh
                      policy.refine(nodes) // nodes of dual mesh --- refined from stencils
                  }
            {
                // define refinement policy
                assert(false);
            }

            auto dual_size() const { return dual_nodes.size(); }

            GridDual(GridDual &&) noexcept = default;
            GridDual(const GridDual &) noexcept = default;
            GridDual() = delete;

            const auto &get_dual_nodes() const
            {
                return dual_nodes;
            }

            const auto operator[](size_t id) const
            {
                return dual_nodes(id);
            }

            // steps between dual nodes
            const DualStepsContainer dual_steps;
            
            // centers of control volumes
            const MeshNodesContainer mesh_nodes; 

        public:
            // dual mesh to be used in simulation
            DualNodesContainer dual_nodes;
            // stencils of the dual mesh.
            // Here, jumps of physical properties occur.
            // These nodes must be included in the dual_mesh_nodes
            // container. So, that operator==() returns true.
            GridDualStencils dual_stencils;

        protected:
            GridDual(
                const GridDualStencils &nodes,
                const DualNodesContainer &refined_mesh) noexcept
                : dual_stencils{nodes}
                , dual_nodes{refined_mesh}
                // make mesh nodes -- centers of control volumes
                , mesh_nodes{CoordinateType::cell_centers(refined_mesh)} 
                , dual_steps{CoordinateType::dual_steps(refined_mesh)
                }
            {
                assert(dual_stencils.size() <= dual_nodes.size());
                auto size{dual_nodes.size()};
                // nodes must be strictly monotonuos
                for (auto idx{size - size}; idx < dual_nodes.size() - 1; ++idx)
                    assert(dual_nodes[idx] < dual_nodes[idx + 1]);
            }
        };

        /// @brief Aggregates grid nodes and
        /// its associated properties, i.e.,
        /// volume per node, heat resistivity etc.
        template <typename CoordinateType_t>
        struct AxesGrid : public GridDual
        {
            using Axes = CoordinateType_t;
        public:
            AxesGrid(const GridDual &dual_nodes) noexcept
                : GridDual{dual_nodes}
                , dual_stencils{dual_nodes.dual_stencils}
                , dual_nodes{dual_nodes.dual_nodes}                                 // copy nodes of dual mesh
                , control_volumes{CoordinateType_t::control_volumes(dual_nodes.dual_nodes)} // make volumes of control cells
                , mesh_steps{CoordinateType_t::mesh_steps(dual_nodes.dual_nodes)}
            {
                assert(dual_nodes.dual_size() > 1ull);
                assert(dual_nodes.dual_size() == dual_steps.size() + 1ull);
                assert(dual_nodes.dual_size() == control_volumes.size() + 1ull);
                assert(dual_nodes.dual_size() == mesh_nodes.size() + 1ull);
                assert(mesh_nodes.size() == mesh_steps.size() + 1ull);
            }

            AxesGrid(AxesGrid &&) noexcept = default;
            AxesGrid(const AxesGrid &) noexcept = default;
            AxesGrid() = delete;

            auto coordinate(auto id) const
            {
                assert(id < mesh_nodes.size());
                return mesh_nodes(id);
            }
            auto dual_coordinate(auto id) const
            {
                assert(id < dual_nodes.size());
                return dual_nodes(id);
            }
            // step between nodes id and id+1
            auto step(auto id) const
            {
                assert(id < mesh_steps.size());
                return mesh_steps(id);
            }

            auto volume(auto idx) const
            {
                assert(idx < control_volumes.size());
                return control_volumes(idx);
            }

            const auto &volumes() const
            {
                return control_volumes;
            }

            auto size() const
            {
                return mesh_nodes.size();
            }

            const GridDualStencils dual_stencils;
            const ControlVolumesContainer control_volumes;
        protected:
            DualNodesContainer dual_nodes;
            // steps between centers of control volumes
            MeshStepsContainer mesh_steps; 
        };

        
        using RGrid = AxesGrid<CoordinateTypes::R_CylCoord>;
        using ZGrid = AxesGrid<CoordinateTypes::Z>;
        
        using XGrid = AxesGrid<CoordinateTypes::X>;
        using YGrid = AxesGrid<CoordinateTypes::Y>;
    } // Grids
} // GPN