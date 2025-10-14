#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>

#include <InjectorDLL/Defines.h>
#include <InjectorDLL/Wrapper.h>

// #include "parse_completion.hpp"

#include <nlohmann/json.hpp>

using VR = std::vector<RealType>;

using namespace std;
using json = nlohmann::json;

int main()
{
    ifstream f("injector_launch.json");
    json data = json::parse(f);

    /*START*/
    // input parameters
    /*fluid*/
    RealType
        viscosity{data["fluid"]["viscosity"].get<RealType>()},
        density{data["fluid"]["density"].get<RealType>()},
        capacity{data["fluid"]["specific_heat_capacity"].get<RealType>()},
        heat_conductivity{data["fluid"]["heat_conductivity"].get<RealType>()},
        joule_thomson{data["fluid"]["joule_thomson"].get<RealType>()};
    /*collector*/
    const VR thickness = data["collector"]["thickness"].get<VR>();
    const VR ext_pressure_stencils{data["collector"]["external_pressure"].get<VR>()};
    // hydrodynamic logs
    const VR is_perforated_stencils = data["collector"]["is_perforated"].get<VR>();
    const VR porosity_stencils = data["collector"]["porosity"].get<VR>();
    const VR permeability_stencils = data["collector"]["permeability"].get<VR>();
    const VR RFP_weights_stencils = data["collector"]["explicit"]["weights"].get<VR>();
    const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
    const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};
    // heat logs
    const VR heatconductivity_stencils = data["collector"]["heatConductivity"].get<VR>();
    const VR solid_density_stencils = data["collector"]["solidDensity"].get<VR>();
    const VR solid_specific_heatcapacity_stencils = data["collector"]["solidSpecificHeatCapacity"].get<VR>();
    /*grid*/
    const RealType
        rMin{data["grid"]["r_start"].get<RealType>()},
        rMax{data["grid"]["r_end"].get<RealType>()},
        q{data["grid"]["r_log_grid"]["q"].get<RealType>()},
        r_max_step{data["grid"]["r_log_grid"]["r_max_step"].get<RealType>()},
        z_minor_step{data["grid"]["z_minor_step"].get<RealType>()}; // m
    /*history*/
    const std::string t_unit = data["history"]["t_unit"].get<std::string>();
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
        t0{factor * data["history"]["start_time"].get<RealType>()};
    RealType t_minor_step {factor * data["history"]["t_minor_step"].get<RealType>()};
    VR t_major_steps = data["history"]["dynamic"]["t_major_step"].get<VR>();
    for (auto &v : t_major_steps)
        v *= factor;
    /*temperatures*/
    const VR well_rates = data["history"]["dynamic"]["well_rate"].get<VR>(); // m^3/s
    const VR inlet_temperatures = data["history"]["dynamic"]["inlet_temperature"].get<VR>();
    /*well*/
    //    const auto casing{parse_completion(data)};
    /*END*/

    const auto &data2 = data["collector"]["geotherma"]["interpolate"];
    const VR geotherma_nodes = data2["z_nodes"].get<VR>();
    const VR geotherma_vals = data2["t_vals"].get<VR>();
    const RealType z_top = data2["z_top"].get<RealType>();

    cout << "Simulation is started." << endl;
    cout << "Please wait..." << endl;

    const auto t_start{chrono::high_resolution_clock::now()};

    Wrapper *instance = new Wrapper(
        // fluid params in SI
        density,           // kg/(m^3)
        capacity,          // J/(kg*K) /* specific heat capacity */
        viscosity,         // Pa*s
        heat_conductivity, // W/(m*K)
        joule_thomson,     // K/bar
        // grid
        rMin,         // m /* typically would be zero */
        rMax,         // m
        q,            // --, q >= 1.0 /* step increment factor */
        r_max_step,   // m /* maximum allowed step in radial direction */
        z_minor_step, // m, /*maximum step within impermeable layers*/
        // eight vectors of the same size
        // values are in SI
        thickness,                 // meter
        ext_pressure_stencils,     // bar
        heatconductivity_stencils, // Watt/(m*K)
        porosity_stencils,         // --
        permeability_stencils,     // m^2
        RFP_weights_stencils,      // -- /*rate distribution between layers of reservoir*/
        is_perforated_stencils,    // {0, 1}, --
        from_coords,               // coordinates of column corrosion, m
        to_layers,
        solid_density_stencils,               // kg/(m^3)
        solid_specific_heatcapacity_stencils, // J/(kg*K)
        // geotherma
        z_top,           // m, /* z-coordinate of the top */
        geotherma_nodes, // m, /* nodes for geotherma interpolation */
        geotherma_vals,  // K, /* reference vals for interpolation */
        // temporal grid
        t0,            // s, start time in seconds
        t_major_steps, // s, in seconds
        t_minor_step,  // s, time step used for numerical integration
        // well
        well_rates,         // ~1.1E-3 m^3/s
        inlet_temperatures, // K
        data);

    const auto t_end{chrono::high_resolution_clock::now()};
    cout << "Elapsed time:                       " << (t_end - t_start).count() * 1E-9 << " seconds\n";

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
