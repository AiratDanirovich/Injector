#pragma once

#include <cmath>

#include <Injector/Grids/Defines.h>

namespace GPN
{
    namespace CoordinateTypes
    {   
        struct GeneralCoordinate
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

        struct CartesianCoordinate : public GeneralCoordinate
        {
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

            /// @brief Interpolate heat conductivity (inverse factor at Laplace term)
            /// @return Heat resistivity at cell face
            static auto face_interpolator(
                RealType xL, RealType xR, 
                RealType xMid,
                RealType valL, RealType valR)
            {
                // coordinates must be monotonous
                assert((xMid-xL)*(xR-xMid) > 0.0);
                assert(xMid != xL);
                assert(xMid != xR);
                assert(xL != xR);

                return (xMid - xL)/valL + (xR - xMid)/valR;
            }
        };

        struct X : public CartesianCoordinate{};
        struct Y : public CartesianCoordinate{};
        struct Z : public CartesianCoordinate{};

        struct R_CylCoord  : public GeneralCoordinate
        {
            static auto control_volumes(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                ControlVolumesContainer out(size);
                
                for(auto idx{size-size}; idx < size; ++idx)
                    out(idx) = volume(nodes(idx), nodes(idx+1));

                return out;
            }
            
            /// @brief Interpolate heat conductivity (inverse factor at Laplace term)
            /// @return Heat resistivity at cell face
            static auto face_interpolator(
                RealType xL, RealType xR, 
                RealType xMid,
                RealType valL, RealType valR)
            {
                // coordinates must be monotonous
                assert((xMid-xL)*(xR-xMid) > 0.0);
                assert(xMid != xL);
                assert(xMid != xR);
                assert(xL != xR);

                return std::log(xMid/xL)/valL + std::log(xR/xMid)/valR;
            }
        protected:
            static RealType volume(RealType x1, RealType x2)
            {
                assert(x2 > x1);
                return (x2 * x2 - x1 * x1) / 2.0;
            }
        };
    } // CoordinateTypes
} // GPN
