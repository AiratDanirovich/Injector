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
    const VR is_perforated_stencils = data["collector"]["is_perforated"];
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
        q{data["grid"]["r_log_grid"]["q"]},
        r_max_step{data["grid"]["r_log_grid"]["r_max_step"]},
        z_minor_step{data["grid"]["z_minor_step"]}; // m
    /*history*/
    const std::string t_unit = data["history"]["t_unit"];
    RealType factor{1.0};
    if (t_unit == "d")
        factor = 24 * 60 * 60;
    else if (t_unit == "h")
        factor = 60 * 60;
    else if (t_unit == "m")
        factor = 60;
    else if (t_unit == "s")
        factor = 1;
    else
        throw std::runtime_error("Incorrect unit of time.");

    const RealType
        t0{data["history"]["start_time"] * factor};
    RealType t_minor_step{data["history"]["t_minor_step"] * factor};
    VR t_major_steps = data["history"]["history_type"]["t_major_step"];
    for (auto &v : t_major_steps)
        v *= factor;
    /*temperatures*/
    const VR well_rates{data["history"]["dynamic"]["well_rate"]}; // m^3/s
    const VR inlet_temperatures{data["history"]["dynamic"]["inlet_temperature"]};
    /*well*/
    const RealType sandface_radius{data["well"]["sandface_radius"]};
    const RealType tube_radius{data["well"]["tube_radius"]};
    /*END*/

    const auto &data2 = data["collector"]["geotherma"]["interpolate"];
    const VR geotherma_nodes = data2["z_nodes"];
    const VR geotherma_vals = data2["t_vals"];
    const RealType z_top = data2["z_top"];

    cout << "Simulation is started." << endl;
    cout << "Please wait..." << endl;

    Wrapper *instance = new Wrapper(
        // fluid params in SI
        density,           // kg/(m^3)
        capacity,          // J/(kg*K) /* specific heat capacity */
        viscosity,         // Pa*s
        heat_conductivity, // W/(m*K)
        // grid
        rMin,         // m /* typically would be zero */
        rMax,         // m
        q,            // --, q >= 1.0 /* step increment factor */
        r_max_step,   // m /* maximum allowed step in radial direction */
        z_minor_step, // m, /*maximum step within impermeable layers*/
        // eight vectors of the same size
        // values are in SI
        thickness,                            // meter
        heatconductivity_stencils,            // Watt/(m*K)
        porosity_stencils,                    // --
        permeability_stencils,                // m^2
        is_permeable_stencils,                // {0, 1}, --
        is_perforated_stencils,               // {0, 1}, --
        solid_density_stencils,               // kg/(m^3)
        solid_specific_heatcapacity_stencils, // J/(kg*K)
        // geotherma
        z_top,           // m, /* z-coordinate of the top */
        geotherma_nodes, // m, /* nodes for geotherma interpolation */
        geotherma_vals,  // K, /* reference vals for interpolation */
        // temporal grid
        t0,           // s, start time in seconds
        t_major_steps,   // s, in seconds
        t_minor_step, // s, time step used for numerical integration
        // well
        tube_radius,      // m
        sandface_radius,  // m
        well_rates,        // ~1.1E-3 m^3/s
        inlet_temperatures // K
    );

    // cout << "After call to DLL\nPress Enter to continue" << endl;
    // getchar();

    const auto &t = instance->get_times();
    cout << "number of saved time moments:       " << t.size() << endl;

    // std::cout << "In main of InjectorUser\nPress Enter to continue" << std::endl;
    // getchar();

    delete instance;

    std::cout << "Simulation is finished.\nPress any key to exit..." << std::endl;
    getchar();

    return 0;
}
