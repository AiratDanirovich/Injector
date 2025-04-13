#pragma once

#include <cmath>

#include <Injector/Grids/Defines.h>

namespace GPN
{
    namespace CoordinateTypes
    {   
        struct StepsCalculator
        {
            static auto steps(const NodesContainer& nodes)
            {
                auto size{nodes.size()};
                assert(size > 0);
                MeshStepsContainer out(size-1);
                for(auto idx{size-size}; idx < size-1; ++idx)
                    out[idx] = nodes[idx+1] - nodes[idx];
                return out;
            }

            // mesh of intermediate, fractional nodes
            static auto dual_mesh(const NodesContainer& nodes)
            {
                auto size{nodes.size()};
                assert(size > 0);
                // dual nodes include domain boundary and intermidiate nodes
                DualNodesContainer out(size+1);

                out(0) = nodes(0);
                for(auto idx{size-size+1}; idx < size+1; ++idx)
                    out(idx) = (nodes(idx-1)+nodes(idx))/2.0;
                auto end{size};
                out(end) = nodes(end-1);

                return out;
            }
            // steps of mesh of intermediate, fractional nodes
            static auto dual_steps(const NodesContainer& nodes)
            {
                auto size{nodes.size()};
                assert(size > 0);
                auto its_dual_mesh{dual_mesh(nodes)};
                // dual nodes include domain boundary and intermidiate nodes
                DualNodesContainer out(size);

                for(auto idx{size-size}; idx < size; ++idx)
                    out(idx) = its_dual_mesh(idx+1)-its_dual_mesh(idx);
                
                return out;
            }
        };

        struct CartesianCoordinate : private StepsCalculator
        {
            static auto volumes(const NodesContainer& nodes)
            {
                auto size{nodes.size()};
                CellVolumeContainer out(size);
                out(0) = (nodes(1) - nodes(0))/2.0;
                for(auto idx{size-size+1}; idx < size-1; ++idx)
                    out(idx) = nodes(idx+1) - nodes(idx);
                out(size-1) = (nodes(size-1) - nodes(size-2))/2.0;

                return out;
            }
            
            using StepsCalculator::steps;
            using StepsCalculator::dual_mesh;
            using StepsCalculator::dual_steps;

        protected:
            // static RealType resistivity(RealType x1, RealType x2)
            // {
            //     return x2 - x1;
            // }
        };

        struct X : CartesianCoordinate{};
        struct Y : CartesianCoordinate{};
        struct Z : CartesianCoordinate{};

        struct RadialCylinderCoordinate  : private StepsCalculator
        {
            static auto volumes(const NodesContainer& nodes)
            {
                auto size{nodes.size()};
                CellVolumeContainer out(size);
                auto mid_val_l{nodes(0)};
                auto mid_val_r{(nodes(1) + nodes(0))/2.0};
                out(0) = volume(nodes(0), mid_val_r);
                
                for(auto idx{size-size+1}; idx < size-1; ++idx)
                {
                    mid_val_l = mid_val_r;
                    mid_val_r = (nodes(idx+1) + nodes(idx))/2.0;
                    out(idx) = volume(mid_val_l, mid_val_r);
                }
                out(size-1) = volume(mid_val_r, nodes(size-1));

                return out;
            }
            
            // steps implementation is inhereted
            using StepsCalculator::steps;
            using StepsCalculator::dual_mesh;
            using StepsCalculator::dual_steps;
        protected:
            static RealType volume(RealType x1, RealType x2)
            {
                assert(x2 > x1);
                return (x2 * x2 - x1 * x1) / 2.0;
            }
            // static RealType resistivity(RealType x1, RealType x2)
            // {
            //     assert(x2 > x1);
            //     assert(x1 > 0.0);
            //     return std::log(x2/x1);
            // }
        };
    } // CoordinateTypes
} // GPN