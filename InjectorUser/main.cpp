#include <iostream>
#include <InjectorDLL/Defines.h>
#include <InjectorDLL/Wrapper.h>

#include <vector>

using VR = std::vector<RealType>;

using namespace std;

int main()
{
    /*START*/
    // input parameters
    /*heat rate*/
    RealType q{1.0}; // not sure about its dimension, just assume = 1.0
    /*fluid*/
    RealType viscosity{6e-4}, density{1000}, capacity{4200};
    /*collector*/
    const RealType rMin{2.0}, rMax{4.0}, zTop{0.0}; // in meters
    const std::ptrdiff_t r_nodes_nmbr{31ull};
    const std::ptrdiff_t nLayers{11ull};
    const VR thickness(nLayers, 0.01); // each layer is 0.01m thick

    const VR conductivity(nLayers, 3.9);
    const VR porosity(nLayers, 1e-16); // porosity assumed zero to exclude contribution of mobile phase, but may be any value < 1.0
    const VR is_permeable(nLayers, 1.0);
    const VR solid_density(nLayers, 2600);
    const VR solid_specific_heatcapacity(nLayers, 770);
    /*temporal grid*/
    const std::ptrdiff_t time_steps_nmbr{51ull};
    const RealType t_start{2e6}; // initial time moment
    const RealType t1{t_start + 2e6};
    const RealType time_step{(t1 - t_start) / time_steps_nmbr};
    const VR time_intervals(time_steps_nmbr, time_step);
    /*END*/

    cout << "before call to DLL\nPress Enter to continue" << endl;
    getchar();

    Wrapper *instance = new Wrapper(
        q, // = 1.0 heat rate
        // fluid params in SI
        density,   // kg/(m^3)
        capacity,  // J/(kg*K) /* specific heat capacity */
        viscosity, // Pa*s
        // grid
        rMin,         // m /* typically would be zero */
        rMax,         // m
        r_nodes_nmbr, // -- /* number of nodes in r-direction, including first and last ones */
        zTop,         // m, /* typically would be zero */
        //    const size_t nZ, // = thickness.size()
        // six vectors of the same size
        // values are in SI
        thickness,                   // meter
        conductivity,                // Watt/(m*K)
        porosity,                    // --
        is_permeable,                // {0, 1}, --
        solid_density,               // kg/(m^3)
        solid_specific_heatcapacity, // J/(kg*K)
        // temporal grid
        t_start, // start time in seconds
                 //    const size_t nt, // = time_intervals.size()
        time_intervals // in seconds
    );

    cout << "After call to DLL\nPress Enter to continue" << endl;
    getchar();

    const auto& t = instance->get_times();

    for(auto i{0ull}; i < t.size(); ++ i)
        std::cout << "t: " << t[i] << std::endl;

    std::cout << "In main of InjectorUser\nPress Enter to exit" << std::endl;
    delete instance;

    getchar();
    return 0;
}
