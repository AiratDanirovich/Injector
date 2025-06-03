#include <iostream>
#include <fstream>
#include <vector>

#include <InjectorDLL/Defines.h>
#include <InjectorDLL/Wrapper.h>

#include <nlohmann/json.hpp>

using VR = std::vector<RealType>;

using namespace std;
using json = nlohmann::json;

VR generate_stencils(RealType t0, RealType t1, RealType t_step_major)
{
    auto segm_count{static_cast<size_t>(std::ceil(t1 - t0) / t_step_major)};
    double step = (t1 - t0) / (segm_count);
    VR out(segm_count + 1ll);

    for (auto i{0ull}; i < out.size(); ++i)
        out[i] = t0 + i * step;
    return out;
}
VR generate_steps(const VR &dual_nodes)
{
  VR out(dual_nodes.size() - 1ll);

  for (auto i{0ull}; i < out.size(); ++i)
    out[i] = dual_nodes[i + 1] - dual_nodes[i];
  return out;
}

int main()
{
    ifstream f("injector_launch.json");
    json data = json::parse(f);

    /*START*/
    // input parameters
    /*fluid*/
    RealType
        viscosity{data["fluid"]["viscosity"]},
        density{data["fluid"]["density"]},
        capacity{data["fluid"]["specificHeatCapacity"]},
        heat_conductivity{data["fluid"]["heatConductivity"]};
    /*collector*/
    const VR thickness = data["collector"]["thickness"];
    //  const ptrdiff_t nLayers{thickness.size()};
    // hydrodynamic logs
    const VR is_permeable_stencils = data["collector"]["is_permeable"];
    const VR porosity_stencils = data["collector"]["porosity"];
    const VR permeability_stencils = data["collector"]["permeability"];
    // heat logs
    const VR heatconductivity_stencils = data["collector"]["heatConductivity"];
    const VR solid_density_stencils = data["collector"]["solidDensity"];
    const VR solid_specific_heatcapacity_stencils = data["collector"]["solidSpecificHeatCapacity"];
    /*grid*/
    const RealType
        rMin{data["grid"]["r_start"]},
        rMax{data["grid"]["r_end"]},
        zTop{data["grid"]["ztop"]}; // m
    const ptrdiff_t rNodes{data["grid"]["rNodes"]};
    /*history*/
    const RealType
        t0{data["history"]["t_start"]},
        t1{data["history"]["t_end"]};
    RealType t_major_step{data["history"]["t_major_step"]};
    RealType t_minor_step{data["history"]["t_minor_step"]};
    t_major_step = std::min(t1 - t0, t_major_step);
    t_minor_step = std::min(t_minor_step, t_major_step);
    const VR t_stencils{generate_stencils(t0, t1, t_major_step)};
    const VR time_steps{generate_steps(t_stencils)};
    /*temperatures*/
    const RealType well_rate{data["history"]["wellRate"]}; // m^3/s
    const RealType initial_temperature{data["collector"]["initTemperature"]};
    const RealType inlet_temperature{data["history"]["inletTemperature"]};
    /*well*/
    const RealType hole_radius{data["well"]["hole_radius"]};
    /*END*/

    // cout << "before call to DLL\nPress Enter to continue" << endl;
    // getchar();

    Wrapper *instance = new Wrapper(
        // fluid params in SI
        density,           // kg/(m^3)
        capacity,          // J/(kg*K) /* specific heat capacity */
        viscosity,         // Pa*s
        heat_conductivity, // W/(m*K)
        // grid
        rMin,   // m /* typically would be zero */
        rMax,   // m
        rNodes, // -- /* number of nodes in r-direction, including first and last ones */
        zTop,   // m, /* typically would be zero */
        // seven +1 vectors of the same size
        // values are in SI
        thickness,                            // meter
        heatconductivity_stencils,            // Watt/(m*K)
        porosity_stencils,                    // --
        permeability_stencils,                // m^2
        is_permeable_stencils,                // {0, 1}, --
        solid_density_stencils,               // kg/(m^3)
        solid_specific_heatcapacity_stencils, // J/(kg*K)
        initial_temperature,                  // K // should be log in the future
        // temporal grid
        t0,           // start time in seconds
                      //    const size_t nt, // = time_intervals.size()
        time_steps,   // in seconds
        t_minor_step, // time step used for numerical integration
        // well
        well_rate,        // ~1.1E-3 m^3/s
        inlet_temperature // K
    );

    // cout << "After call to DLL\nPress Enter to continue" << endl;
    // getchar();

    const auto &t = instance->get_times();
    cout << "number of saved time moments:       " << t.size() << endl;

    // std::cout << "In main of InjectorUser\nPress Enter to continue" << std::endl;
    // getchar();

    delete instance;

    std::cout << "Simulation done\nPress Enter to exit" << std::endl;
    //getchar();

    return 0;
}
