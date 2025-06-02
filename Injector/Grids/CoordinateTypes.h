#pragma once

#include <cmath>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/CoordinateSystem.hpp>

namespace GPN
{
    namespace CoordinateTypes
    {   
        /// @brief Calculations associated with any coordinate axes,
        /// i.e., steps between (dual and regular) adjacent nodes, 
        /// center coordinates between two nodes
        struct GeneralCoordinate
        {
            /// @brief Normal distance between two faces of control volume
            static auto dual_steps(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                assert(size > 0);
                DualStepsContainer out(size);
                for(auto idx{size-size}; idx < size; ++idx)
                    out(idx) = nodes(idx+1) - nodes(idx);
                return out;
            }

            /// @brief 
            /// @param nodes Nodes of dual mesh
            /// @return Centers of control volumes
            static auto cell_centers(const DualNodesContainer& nodes)
            {
                assert(nodes.size() > 2ll);

                auto size{nodes.size()-1ll};
                MeshNodesContainer out(size);
                // centers of boundary cells are moved to the domain boundary
                for(auto idx{0ll}; idx < size-0ll; ++idx)
                    out(idx) = (nodes(idx+1) + nodes(idx))/2.0;

                const RealType tol = 1e-12;
                out.head(1ll) = nodes.head(1ll)+tol;
                out.tail(1ll) = nodes.tail(1ll)-tol;

                return out;
            }

            /// @brief 
            /// @param nodes Nodes of dual mesh
            /// @return Steps between centers of control volumes
            static auto mesh_steps(const DualNodesContainer& nodes)
            {
                auto mesh_nodes{cell_centers(nodes)};

                assert(mesh_nodes.size() > 1ll);
                auto size{mesh_nodes.size()-1ll};
                MeshStepsContainer out(size);
                for(auto id{0ll}; id < size; ++id)
                    out(id) = mesh_nodes(id+1)-mesh_nodes(id);
                return out;                    
            }
        };

        /// @brief Calculations associated with Cartesian coordinate
        struct CartesianCoordinate : public GeneralCoordinate
        {
            /// @brief Generate control volumes from dual mesh
            /// @param nodes Nodes of dual mesh
            /// @return Volumes of control cells
            static auto control_volumes(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                ControlVolumesContainer out(size);
                for(auto idx{0ll}; idx < (ptrdiff_t)size; ++idx)
                    out(idx) = nodes(idx+1) - nodes(idx);
                return out;
            }

            /// @brief Interpolate heat conductivity (factor at Laplace term)
            /// @return Heat conductivity at cell face
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

                if((valL == 0.0) || (valR == 0.0))
                    return 0.0;

                return 1.0/((xMid - xL)/valL + (xR - xMid)/valR);
            }

            
            /// @brief Interpolate const heat conductivity (factor at Laplace term)
            /// @return Heat conductivity at cell face
            static auto const_face_interpolator(
                RealType xL,
                RealType xR,
                const Eigen::ArrayX<RealType>& val)
            {
                assert(xL != xR);
                return 1.0/((xR-xL)/val);
            }
        };

        /// @brief X coordinate
        struct X : public CartesianCoordinate{};
        /// @brief Y coordinate
        struct Y : public CartesianCoordinate{};
        /// @brief Z coordinate
        struct Z : public CartesianCoordinate{};
        
        /// @brief R coordinate of cylinder (polar) system of coordinates
        struct R_CylCoord  : public GeneralCoordinate
        {
            /// @brief Calculates volumes of grid cells for radial coordinate
            /// @param nodes Nodes of dual mesh
            /// @return Volumes of control volume cells
            static auto control_volumes(const DualNodesContainer& nodes)
            {
                auto size{nodes.size()-1ull};
                ControlVolumesContainer out(size);
                
                for(auto idx{size-size}; idx < size; ++idx)
                    out(idx) = volume(nodes(idx), nodes(idx+1));

                return out;
            }
            
            /// @brief Interpolate heat conductivity (factor at Laplace term)
            /// @return Heat conductivity at cell face
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

                if((valL == 0.0) || (valR == 0.0))
                    return 0.0;
                    
                return 1.0/(std::log(xMid/xL)/valL + std::log(xR/xMid)/valR);
            }
            /// @brief Interpolate const heat conductivity (factor at Laplace term)
            /// @return Heat conductivity at cell face
            static auto const_face_interpolator(
                RealType xL,
                RealType xR,
                const Eigen::ArrayX<RealType>& val)
            {
                assert(xL != xR);
                return 1.0/(std::log(xR/xL)/val);
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
