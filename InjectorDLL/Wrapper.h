#pragma once

#include <vector>

#ifdef MYLIBRARY_EXPORT
#define LIBRARY_API __declspec(dllexport)
#else
#define LIBRARY_API __declspec(dllimport)
#endif

#include <vector>
#include <InjectorDLL/Defines.h>

using VR = std::vector<RealType>;

class CustomVector;
#pragma warning(push)
#pragma warning(disable : 4251)
class LIBRARY_API Wrapper
{
public:
    Wrapper(
        // fluid params in SI
        const RealType density,                 // kg/(m^3)
        const RealType capacity,                // J/(kg*K) /* specific heat capacity */
        const RealType viscosity,               // Pa*s
        const RealType heat_conductivity_fluid, // W/(m*K)
        // grid
        const RealType rMin,         // m /* typically would be zero */
        const RealType rMax,         // m
        const RealType q,            // --, q >= 1.0 /* step increment factor */
        const RealType r_max_step,   // m /* maximum allowed step in radial direction */
        const RealType zTop,         // m, /* typically would be zero */
        const RealType z_minor_step, // m, /*maximum step within impermeable layers*/
                                     // eight +1 vectors of the same size
        // values are in SI
        const VR &thickness,                   // meter
        const VR &conductivity,                // W/(m*K)
        const VR &porosity,                    // --
        const VR &permeability,                // m^2
        const VR &is_permeable,                // {0, 1}, --
        const VR &is_perforated,               // {0, 1}, --
        const VR &solid_density,               // kg/(m^3)
        const VR &solid_specific_heatcapacity, // J/(kg*K)
        const RealType initial_temperature,    // K // should be log in the future
        // temporal grid
        const RealType t_start,      // start time in seconds
        const VR &time_intervals,    // time intervals (in seconds) of const rates
        const RealType t_minor_step, // time step used for numerical integration
        // well
        const RealType tube_radius,      // m
        const RealType sandface_radius,  // m
        const RealType well_rate,        // ~1.1E-3 m^3/s
        const RealType inlet_temperature // K
    );
    ~Wrapper();

    std::vector<RealType> get_times() const;
    // std::vector<std::vector<RealType>> get_temps() const;

private:
    std::vector<std::vector<RealType>> t_radial_distribution;
    std::vector<RealType> time;
};
#pragma warning(pop)