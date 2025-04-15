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

            auto size() const
            {
                return dual_nodes.size();
            }

        protected:
            // read-only
            DualNodesContainer dual_nodes;
        };

        /// @brief Container for dual grid nodes (refined as well as stencils).
        /// The nodes are somehow distributed, uniformly or non-uniformly, between adjuscent stencils.
        struct GridDual
        {
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

            auto size() const { return dual_nodes.size(); }

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
                : dual_stencils{nodes}, dual_nodes{refined_mesh}
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
        struct Grid1D
        {
        public:
            Grid1D(const GridDual &dual_nodes) noexcept
                : dual_nodes{dual_nodes.dual_nodes}                                 // copy nodes of dual mesh
                , dual_steps{CoordinateType_t::dual_steps(dual_nodes.dual_nodes)} // normal distance between two faces of control volume
                , control_volumes{CoordinateType_t::control_volumes(dual_nodes.dual_nodes)} // make volumes of control cells
                , mesh_nodes{CoordinateType_t::cell_centers(dual_nodes.dual_nodes)} // make mesh nodes -- centers of control volumes
                , mesh_steps{CoordinateType_t::mesh_steps(dual_nodes.dual_nodes)}
            {
                assert(dual_nodes.size() > 1ull);
                assert(dual_nodes.size() == dual_steps.size() + 1ull);
                assert(dual_nodes.size() == control_volumes.size() + 1ull);
                assert(dual_nodes.size() == mesh_nodes.size() + 1ull);
                assert(mesh_nodes.size() == mesh_steps.size() + 1ull);
            }

            Grid1D(Grid1D &&) noexcept = default;
            Grid1D(const Grid1D &) noexcept = default;
            Grid1D() = delete;

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

            auto dual_size() const
            {
                return dual_nodes.size();
            }

        protected:
            DualNodesContainer dual_nodes;
            DualStepsContainer dual_steps;
            ControlVolumesContainer control_volumes;
            MeshNodesContainer mesh_nodes; // centers of control volumes
            MeshStepsContainer mesh_steps; // steps between centers of control volumes
        };

        /// @brief Two-dimensional grid
        /// @tparam FirstDir Type for grid in first [r] direction
        /// @tparam SecondDir Type for grid in second [z] direction
        template <typename FirstDir, typename SecondDir>
        struct StructuredGrid2D
        {
        public:
            struct Point
            {
                RealType x, y;
            };

            FirstDir first_coord;
            SecondDir second_coord;

            StructuredGrid2D(
                const FirstDir &first_coord, 
                const SecondDir &second_coord)
                    : first_coord{first_coord}
                    , second_coord{second_coord}
                    , its_volumes(
                        first_coord.size(),
                        second_coord.size()
                    )
            {
                // set volumes
                for (std::ptrdiff_t j = 0; j < its_volumes.cols(); ++j)
                    for (std::ptrdiff_t i = 0; i < its_volumes.rows(); ++i)
                        its_volumes(i, j) =
                            first_coord.volume(i) *
                            second_coord.volume(j);
            }

            auto ccordinates(auto id1, auto id2) const
            {
                return Point{first_coord.mesh_nodes(id1), second_coord.mesh_nodes(id2)};
            }

            // steps in two directions,
            // between nodes id1 and id1+1 in first direction
            // and between nodes id2 and id2+1 in second direction
            auto steps(auto id1, auto id2) const
            {
                return {first_coord.mesh_steps(id1), second_coord.mesh_steps(id2)};
            }

            /// @brief cell volume at node ids {id1, id2}
            auto volume(auto id1, auto id2) const
            {
                return first_coord.cell_volumes(id1) * second_coord.cell_volumes(id2);
            }

            const auto &volumes() const
            {
                return its_volumes;
            }

            CellVolumeContainer2D its_volumes;
        };

        struct StructuredCylinderGrid2D
            : public StructuredGrid2D<
                  Grid1D<CoordinateTypes::Z>,
                  Grid1D<CoordinateTypes::RadialCylinderCoordinate>>
        {
            using StructuredGrid2D<
                Grid1D<CoordinateTypes::Z>,
                Grid1D<CoordinateTypes::RadialCylinderCoordinate>>::StructuredGrid2D;
        };

    } // Grids
} // GPN