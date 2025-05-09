#pragma once

#include <vector>

#ifdef MYLIBRARY_EXPORT
#define LIBRARY_API __declspec(dllexport)
#else
#define LIBRARY_API __declspec(dllimport)
#endif

#include <vector>
#include <Defines.h>

using VR = std::vector<RealType>;

class CustomVector;
#pragma warning( push )
#pragma warning( disable : 4251 )
class LIBRARY_API Wrapper
{
public:
    Wrapper(
        const RealType q, // = 1.0 heat rate
        // fluid params in SI
        const RealType density,   // kg/(m^3)
        const RealType capacity,  // J/(kg*K) /* specific heat capacity */
        const RealType viscosity, // Pa*s
        // grid
        const RealType rMin,       // m /* typically would be zero */
        const RealType rMax,       // m
        const size_t r_nodes_nmbr, // -- /* number of nodes in r-direction, including first and last ones */
        const RealType zTop,       // m, /* typically would be zero */
        //    const size_t nZ, // = thickness.size()
        // six vectors of the same size
        // values are in SI
        const VR &thickness,                   // meter
        const VR &conductivity,                // Watt/(m*K)
        const VR &porosity,                    // --
        const VR &is_permeable,                // {0, 1}, --
        const VR &solid_density,               // kg/(m^3)
        const VR &solid_specific_heatcapacity, // J/(kg*K)
        // temporal grid
        const RealType t_start,  // start time in seconds
                                 //    const size_t nt, // = time_intervals.size()
        const VR &time_intervals // time intervals (in seconds) of const rates
    );
    ~Wrapper();

    std::vector<RealType> get_times() const;
    std::vector<std::vector<RealType>> get_temps() const;

private:
    std::vector<std::vector<RealType>> t_radial_distribution;
    std::vector<RealType> time;
};
#pragma warning( pop )