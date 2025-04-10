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
                    assert(nodes[idx] > nodes[idx + 1]);

                // copy stencils to local container
                for (size_t idx{0}; idx < nodes.size(); ++idx)
                    mesh_nodes[idx] = nodes[idx];
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
            // nodes must be strictly monotonuos
            for (std::ptrdiff_t idx{0}; idx < mesh_nodes.size() - 1; ++idx)
                assert(mesh_nodes[idx] > mesh_nodes[idx + 1]);
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
            , mesh_steps{CoordinateType_t::steps(nodes)}
            , cell_volumes{CoordinateType_t::volumes(nodes)}
            {
                assert(nodes.size() > 1ull);
            }

            Grid1D(Grid1D &&) noexcept = default;
            Grid1D(const Grid1D &) noexcept = default;
            Grid1D() = delete;

            // step between nodes id and id+1
            auto step(size_t id) const
            {
                assert(id < mesh_steps.size());
                return mesh_steps(id);
            }

            auto cell_volume(size_t idx) const{
                assert(idx < cell_volumes.size());
                return cell_volumes(idx);
            }

        protected:
            GridNodes mesh_nodes;
            MeshStepsContainer mesh_steps;
            CellVolumeContainer cell_volumes;
        };

        template<typename CoordinateType_t>
        struct Grid1DUniform : public Grid1D<CoordinateType_t>
        {
        };

        
        /// @brief To generate cell volume along the
        /// Cartsian direction.
        template<typename CoordinateType_t>
        struct GridCartesianCoord : public Grid1D<CoordinateType_t>
        {
        };


        /// @brief To generate cell volume along the
        /// radial direction of cylinrical coordinates
        template<typename CoordinateType_t>
        struct GridRadialCoord : public Grid1DUniform<CoordinateType_t>
        {
        };

        struct UniformGrid1D : private std::vector<float_t>
        {
            using Cell_Volume = Eigen::ArrayX<float_t>;

        public:
            static UniformGrid1D CreateFromStep(float_t a, float_t b, float_t step)
            {
                assert(b > a);
                assert(step > 0.0);
                assert(step < b - a);

                size_t segm_nmbr = static_cast<size_t>((b - a) / step) + 1;

                return CreateFromNodes(a, b, segm_nmbr + 1);
            }

            static UniformGrid1D CreateFromNodes(float_t a, float_t b, size_t nodes_nmbr)
            {
                assert(b > a);
                assert(nodes_nmbr > 1ull);

                float_t step = (b - a) / (nodes_nmbr - 1ull);
                std::vector<float_t> nodes(nodes_nmbr, a);

                std::generate(nodes.begin(), nodes.end(), [n = 0, &step, &a]() mutable
                              { return a + n++ * step; });

                return {nodes};
            }

            UniformGrid1D(const UniformGrid1D &) noexcept = default;
            UniformGrid1D(UniformGrid1D &&) noexcept = default;

            float_t step(ptrdiff_t) const
            {
                return its_step;
            }

            const Cell_Volume &cell_volume() const
            {
                return its_cell_volume;
            }

        public:
            using std::vector<float_t>::operator[];
            using std::vector<float_t>::data;
            using std::vector<float_t>::size;
            using std::vector<float_t>::begin;
            using std::vector<float_t>::end;
            using std::vector<float_t>::front;
            using std::vector<float_t>::back;
            using std::vector<float_t>::vector;

        protected:
            UniformGrid1D(const std::vector<float_t> &nodes)
                : std::vector<float_t>{nodes}
            {
                assert(nodes.size() > 1);
                its_step = nodes[1] - nodes[0];

                // volumes of cells around every node
                // a Uniform grid is assumed
                its_cell_volume = Cell_Volume::Constant(nodes.size(), its_step);
                its_cell_volume(0) /= 2.0;
                its_cell_volume(its_cell_volume.size() - 1) /= 2.0;
            }

            float_t its_step;

            Cell_Volume its_cell_volume;
        };
    } // Grids
} // GPN