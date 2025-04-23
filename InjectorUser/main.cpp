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
    RealType q{1.0};
    /*fluid*/
    RealType viscosity{6e-4}, density{1000}, capacity{4200};
    /*collector*/
    const RealType rMin{1.0}, rMax{2.0}, zTop{0.0};
    const std::ptrdiff_t r_nodes_nmbr{31ull};
    const std::ptrdiff_t nLayers{11ull};
    const VR thickness(nLayers, 0.01); // each layer is 1m thick

    const VR conductivity(nLayers, 3.9);
    const VR porosity(nLayers, 1e-16);
    const VR is_permeable(nLayers, 1.0);
    const VR solid_density(nLayers, 3.9 /*should be 2600 in SI*/);
    const VR solid_specific_heatcapacity(nLayers, 1.0 /*should be 770 in SI*/);
    /*temporal grid*/
    const std::ptrdiff_t time_steps_nmbr{51ull};
    const RealType t_start{1.0}; // initial time moment
    const RealType t1{t_start + 1.0};
    const RealType time_step{(t1 - t_start) / time_steps_nmbr};
    const VR time_intervals(time_steps_nmbr, time_step);
    /*END*/

    cout << "before call to DLL" << endl;
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
        time_intervals);

    cout << "After call to DLL" << endl;
    getchar();

    const auto& t = instance->get_times();

    for(auto i{0ull}; i < t.size(); ++ i)
        std::cout << "t: " << t[i] << std::endl;

    std::cout << "In main of InjectorUser" << std::endl;
    delete instance;

    getchar();
    return 0;
}
