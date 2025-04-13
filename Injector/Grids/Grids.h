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
        /// @brief Struct defines the grid stencils, i.e.
        /// nodes that must be included as grid nodes.
        /// Further refinement is possible
        struct GridStencils
        {
            GridStencils(const std::vector<RealType> &nodes) noexcept 
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
            // Here, jumps of physical properties occur.
            // These nodes must be included in th mesh_nodes
            // container. So, that operator== returned true.
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
        public:
            Grid1D(const GridNodes &nodes) noexcept 
            : mesh_nodes{nodes}
            , mesh_steps{CoordinateType_t::steps(nodes.mesh_nodes)}
            , cell_volumes{CoordinateType_t::volumes(nodes.mesh_nodes)}
            , dual_nodes{CoordinateType_t::dual_nodes(nodes.mesh_nodes)}
            , dual_steps{CoordinateType_t::dual_steps(nodes.mesh_nodes)}
            {
                assert(nodes.size() > 1ull);
            }

            Grid1D(Grid1D &&) noexcept = default;
            Grid1D(const Grid1D &) noexcept = default;
            Grid1D() = delete;

            auto coord(auto id) const
            {
                assert(id < mesh_nodes.size());
                return mesh_nodes.get_nodes()(id);
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

            const auto& volumes() const
            {return cell_volumes;}

            auto size() const
            {return mesh_nodes.size();}

        protected:
            GridNodes mesh_nodes;
            // distance between adjacent nodes
            MeshStepsContainer mesh_steps;
            CellVolumeContainer cell_volumes;

            DualStepsContainer dual_steps;
            DualNodesContainer dual_nodes;
        };

        /// @brief Two-dimensional grid
        /// @tparam FirstDir Type for grid in first [r] direction
        /// @tparam SecondDir Type for grid in second [z] direction
        template<typename FirstDir, typename SecondDir>
        struct StructuredGrid2D
        {
        public:
            struct Point
            {
                RealType x,y;
            };

            FirstDir first_coord;
            SecondDir second_coord;

            StructuredGrid2D(const FirstDir& first_coord, const SecondDir& second_coord) 
                : first_coord{first_coord}
                , second_coord{second_coord}
                , its_volumes(
                    first_coord.size(),
                    second_coord.size())
            {
                // set volumes
                for (std::ptrdiff_t j = 0; j < its_volumes.cols(); ++j)
                    for (std::ptrdiff_t i = 0; i < its_volumes.rows(); ++i)
                        its_volumes(i, j) =
                                first_coord.volume(i)*
                                second_coord.volume(j);
            }

            auto coords(auto id1, auto id2) const
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
            auto volume(auto id1, auto id2) const{
                return first_coord.cell_volumes(id1)*second_coord.cell_volumes(id2);
            }

            const auto& volumes() const{
                return its_volumes;
            }

            CellVolumeContainer2D its_volumes;
        };

        struct StructuredCylinderGrid2D 
        : public StructuredGrid2D<Grid1D<CoordinateTypes::Z>, Grid1D<CoordinateTypes::RadialCylinderCoordinate>>
        {
            using StructuredGrid2D<Grid1D<CoordinateTypes::Z>, Grid1D<CoordinateTypes::RadialCylinderCoordinate>>::StructuredGrid2D;
        };

    } // Grids
} // GPN