#pragma once

#include <cmath>

#include <Injector/Grids/Defines.h>

namespace GPN
{
    namespace CoordinateTypes
    {   
        struct CartesianCoordinate
        {
            /// @brief Normal distance between two faces of control volume
            static auto dual_steps(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                assert(size > 0);
                DualStepsContainer out(size);
                for(auto idx{size-size}; idx < size-1; ++idx)
                    out(idx) = nodes(idx+1) - nodes(idx);
                return out;
            }
            
            /// @brief Generate control volumes from dual mesh
            /// @param nodes 
            /// @return 
            static auto control_volumes(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                ControlVolumesContainer out(size);
                for(auto idx{size-size}; idx < size; ++idx)
                    out(idx) = nodes(idx+1) - nodes(idx);
                return out;
            }

            static auto cell_centers(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                MeshNodesContainer out(size);
                for(auto idx{size-size}; idx < size; ++idx)
                    out(idx) = (nodes(idx+1) + nodes(idx))/2.0;
                return out;
            }

            static auto mesh_steps(const DualNodesContainer& nodes)
            {
                auto mesh_nodes{cell_centers(nodes)};

                assert(mesh_nodes.size() >= 1);
                auto size{mesh_nodes.size()-1ull};
                MeshStepsContainer out(size);
                for(auto id{size-size}; id < size; ++id)
                    out(id) = mesh_nodes(id+1)-mesh_nodes(id);
                return out;                    
            }
        };

        struct X : public CartesianCoordinate{};
        struct Y : public CartesianCoordinate{};
        struct Z : public CartesianCoordinate{};

        struct RadialCylinderCoordinate  : private CartesianCoordinate
        {
            using CartesianCoordinate::dual_steps;

            static auto control_volumes(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                ControlVolumesContainer out(size);
                
                for(auto idx{size-size}; idx < size; ++idx)
                    out(idx) = volume(nodes(idx), nodes(idx+1));

                return out;
            }
            
            using CartesianCoordinate::cell_centers;
            using CartesianCoordinate::mesh_steps;
        protected:
            static RealType volume(RealType x1, RealType x2)
            {
                assert(x2 > x1);
                return (x2 * x2 - x1 * x1) / 2.0;
            }
        };
    } // CoordinateTypes
} // GPN