#include <vector>
#include <cmath>
#include <memory>
#include <numbers>
#include <fstream>
#include <iostream>
#include <filesystem>

#include <InjectorDLL/Wrapper.h>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/Grids2D.hpp>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/History/History.hpp>
#include <Injector/History/RatesFactory.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Collector.hpp>
#include <Injector/Model/WellFactory.hpp>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Solver/BoundaryConditions.hpp>
#include <Injector/Solver/State2D.hpp>
#include <Injector/Solver/InitialCondition.hpp>
#include <Injector/Solver/FullImplicit/Solver.hpp>
#include <Injector/Solver/SolverManager.hpp>

#include <Eigen/Core>

#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

using namespace GPN;
using namespace GPN::Phases;
using namespace GPN::Completion;
using namespace GPN::EqSolver;
using namespace GPN::EqSolver::FullImplicit;

namespace fs = std::filesystem;

/// @brief Initial temperature is assumed to be constant
struct FunctorIC : public InitialConditions::ICFunctorBase
{
    FunctorIC(const Logs::Geotherma &geotherma)
        : geotherma{geotherma}
    {
    }
    RealType operator()(const ptrdiff_t z_id, const ptrdiff_t, RealType) const override
    {
        return geotherma(z_id);
    }

protected:
    const Logs::Geotherma &geotherma;
};

template <typename Grid_t_ptr>
auto ICFactory(RealType t_start, const Grid_t_ptr grid, const Logs::Geotherma &geotherma)
{
    return State::State2D{State::State2D::FillWithFunctor(*grid, FunctorIC{geotherma}, t_start)};
}

template <typename Well_t>
struct FunctorBC : public GPN::BoundaryConditions::BCFunctorBase
{
    using Grid2D_t = Grids::StructuredCylinderGrid2DAxisymmetric;
    using ConvectionFieldFactory_t =
        GPN::FaceProperties::RatesFactory<
            Grid2D_t, Well_t, PhaseProperties>;

    FunctorBC(
        const Logs::IsPermeable &is_permeable,
        const ConvectionFieldFactory_t &flow_field, // volumetric flow rate
        const cptr<const Grid2D_t> grid_ptr)
        : flow_field{flow_field},
          is_permeable{is_permeable},
          grid_ptr{grid_ptr}
    {
    }

    RealType operator()(ptrdiff_t z_id, RealType r, RealType t) const override
    {
        if (r == grid_ptr->second_coord.dual_front())
            return flow_field.get_flow_in_axes2()(z_id, 0ll) * flow_field.get_temperature();

        if (r == grid_ptr->second_coord.dual_back())
            return 0.0;

        assert(false);
        return 0.0;
    }

    RealType operator()(RealType z, ptrdiff_t r_id, RealType t) const override
    {
        if (z == grid_ptr->first_coord.dual_front())
            return flow_field.get_flow_in_axes1()(0ll, r_id) * flow_field.get_temperature();

        if (z == grid_ptr->first_coord.dual_back())
            return 0.0;

        assert(false);
        return 0.0;
    }

protected:
    const ConvectionFieldFactory_t &flow_field;
    const Logs::IsPermeable &is_permeable;
    const cptr<const Grid2D_t> grid_ptr;
};

const auto make_completion(const std::array<std::array<RealType, 6>, 6> &data)
{
    using namespace GPN::Completion;

    std::vector<Ring> out;
    out.reserve(6);

    { // flowing fluid
        const auto &data2 = data[0ull];
        out.push_back(
            Ring{
                Flow{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{data2[5ull]}});
    }

    { // tube
        const auto &data2 = data[1ull];
        out.push_back(
            Ring{
                Tube{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{data2[5ull]}});
    }

    { // annulus
        const auto &data2 = data[2ull];
        out.push_back(
            Ring{
                Annulus{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{data2[5ull]}});
    }

    { // column
        const auto &data2 = data[3ull];
        out.push_back(
            Ring{
                Column{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{data2[5ull]}});
    }

    { // cement_1
        const auto &data2 = data[4ull];
        out.push_back(
            Ring{
                Cement{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{data2[5ull]}});
    }

    { // cement_2
        const auto &data2 = data[5ull];
        out.push_back(
            Ring{
                Cement{
                    Density{data2[0ull]},
                    SpecificHeatCapacity{data2[1ull]},
                    GPN::HeatConductivity{data2[2ull]}},
                Thickness{data2[3ull]},
                InnerRadius{data2[4ull]},
                Depth{data2[5ull]}});
    }

    return out;
}

Wrapper::Wrapper(
    // fluid params in SI
    const RealType density,                 // kg/(m^3)
    const RealType capacity,                // J/(kg*K) /* specific heat capacity */
    const RealType viscosity,               // Pa*s
    const RealType heat_conductivity_fluid, // Watt/(m*K)
    // grid
    const RealType rMin,         // m /* typically would be zero */
    const RealType rMax,         // m
    const RealType q,            // --, q >= 1.0 /* step increment factor */
    const RealType r_max_step,   // m /* maximum allowed step in radial direction */
    const RealType z_minor_step, // m, /*maximum step within impermeable layers*/
    // eight vectors of the same size
    // values are in SI
    const VR &thickness,                   // meter
    const VR &solid_heatconductivity,      // Watt/(m*K)
    const VR &porosity,                    // 0.0 < porosity <= 1.0, --
    const VR &permeability_stencils,       // m^2
    const VR &weights_stencils,            // -- /*rate distribution between layers*/
    const VR &is_permeable,                // {0, 1}, --
    const VR &is_perforated,               // {0, 1}, --
    const VR &solid_density,               // kg/(m^3)
    const VR &solid_specific_heatcapacity, // J/(kg*K)
    // geotherma
    const RealType z_top,      // m, /* z-coordinate of the top */
    const VR &geotherma_nodes, // m, /* nodes for geotherma interpolation */
    const VR &geotherma_vals,  // K, /* reference vals for interpolation */
    // temporal grid
    const RealType t_start,      // start time in seconds
    const VR &time_intervals,    // intervals of const rates)
    const RealType t_minor_step, // time step used for numerical integration
    // well
    const VR &well_rates,         // ~1.1E-3 m^3/s
    const VR &inlet_temperatures, // K
    // casing
    // {fluid, tube, annulus, column, cementInner, cementOuter}
    const std::array<MaterialProps, 6> &casing_data)
{
    // adapt stl container to Eigen container
    LogValuesContainer is_permeable_stencils(is_permeable.size());
    std::copy(is_permeable.cbegin(), is_permeable.cend(), is_permeable_stencils.begin());
    LogValuesContainer is_perforated_stencils(is_perforated.size());
    std::copy(is_perforated.cbegin(), is_perforated.cend(), is_perforated_stencils.begin());
    LogValuesContainer solid_density_stencils(solid_density.size());
    std::copy(solid_density.begin(), solid_density.end(), solid_density_stencils.begin());
    LogValuesContainer solid_specific_heatcapacity_stencils(solid_specific_heatcapacity.size());
    std::copy(solid_specific_heatcapacity.begin(), solid_specific_heatcapacity.end(), solid_specific_heatcapacity_stencils.begin());
    LogValuesContainer porosity_stencils(porosity.size());
    std::copy(porosity.begin(), porosity.end(), porosity_stencils.begin());
    LogValuesContainer solid_heatconductivity_stencils(solid_heatconductivity.size());
    std::copy(solid_heatconductivity.begin(), solid_heatconductivity.end(), solid_heatconductivity_stencils.begin());

    const Casing completion{make_completion(casing_data)};

    // r_stencils
    GPN::WellHoles well_holes{WellHolesFactory::create(completion)};
    VR r_stencils = well_holes.generate_log_radial_grid(
        rMin, rMax, q, r_max_step);
    // z-refiner
    GPN::Grids::RefinerVerticle refiner{z_minor_step, is_permeable_stencils};
    // the grid itself
    const auto grid2D{
        Grids::CylinderGridFactory::create(
            refiner,
            Grids::Factory::generate_dual_grid_stencils_from_steps(
                0.0, thickness),
            r_stencils)};
    const auto &grid_r{grid2D->second_coord};
    const auto &grid_z{grid2D->first_coord};
    // collector
    const Logs::Rocks::CoreSampleLogs core_data{
        is_permeable_stencils,
        is_perforated_stencils,
        porosity_stencils,
        permeability_stencils,
        grid_z};
    // make fluid
    const PhaseProperties water{
        FluidFactory::create_water(
            Viscosity{viscosity},
            Density{density},
            SpecificHeatCapacity{capacity},
            HeatConductivity{heat_conductivity_fluid})};

    const auto weights{
        Logs::RateWeightsFactory::create(
            weights_stencils,
            core_data.is_permeable,
            grid_z)};

    const Well_Explicit well{
        core_data.is_permeable,
        core_data.is_perforated, weights};

    const Logs::Rocks::HeatLogs heat_logs{
        solid_density_stencils,
        solid_specific_heatcapacity_stencils,
        solid_heatconductivity_stencils,
        porosity_stencils,
        water,
        grid2D->first_coord};

    Properties::Rocks::HeatProps heat_props{
        heat_logs, grid2D};
    heat_props.apply_well(completion, well);

    FaceProperties::Rocks::HeatFaceProps heat_face_props{
        heat_props, grid2D};
    heat_face_props.apply_well(completion, well);

    // history
    const History history{
        HistoryFactory::createFixedRate(time_intervals, well_rates, inlet_temperatures)};
    // rates field factory
    FaceProperties::RatesFactory rates_factory{
        grid2D, well, history, water};
    // initial condition
    Logs::Geotherma geotherma{
        Logs::GeothermaFactory::create(
            geotherma_nodes,
            geotherma_vals,
            z_top,
            grid2D->first_coord)};
    const auto initial_state{ICFactory(t_start, grid2D, geotherma)};
    // boundary conditions
    const GPN::BoundaryConditions::BoundaryConditions bc{
        *grid2D,
        std::make_shared<FunctorBC<std::remove_const<decltype(well)>::type>>(
            core_data.is_permeable, rates_factory, grid2D),
        BoundaryConditions::BoundaryCondition::second};
    // solver
    using Solver_t = decltype(Solver{
        heat_face_props.medium_heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        rates_factory, initial_state,
        bc, t_start});

    auto solver_ptr = std::make_shared<Solver_t>(
        heat_face_props.medium_heat_conductivity,
        grid2D,
        heat_props.medium_vol_heatcapacity,
        rates_factory, initial_state,
        bc, t_start);

    const auto &solver{*solver_ptr};

    SolverManager solver_manager{history, solver_ptr};

    solver_manager.run(t_minor_step);

    const auto &[times, states] = solver.solution();
    this->time = times;

    if (!fs::is_directory("output") || !fs::exists("output")) // Check if src folder exists
    {
        fs::create_directory("output"); // create src folder
    }

    try
    {
        ifstream f("separators.json");
        json data = json::parse(f);

        const std::string sep = data["coeff_sep"];

        auto grid{grid2D->second_coord.mesh_nodes};
        grid(0ll) = grid_r.dual_nodes(1ll);

        const Eigen::IOFormat commaFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, sep, sep, "", "", "", "");
        for (auto z{0ll}, layer_id{0ll}; z < grid2D->first_coord.mesh_nodes.size(); ++z)
        {
            if (core_data.is_permeable(z) == 1.0)
            {
                ofstream f{std::string{"output/layer_"} + std::to_string(layer_id) + std::string{".csv"}};

                f << sep << sep << grid.transpose().format(commaFmt) << '\n';
                for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
                {
                    f << t << sep << times[t] << sep << states[t].cur_state.row(z).format(commaFmt) << '\n';
                }

                f.close();
                ++layer_id;
            }
        }

        {
            ofstream f{std::string{"output/well_temperature.csv"}};
            f << sep << sep << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                f << t << sep << times[t] << sep << states[t].cur_state.col(0ll).format(commaFmt) << '\n';
            }
            f.close();
        }

        {
            ofstream f{std::string{"output/cement_temperature.csv"}};
            f << sep << sep << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            for (auto t{0ll}; t < (ptrdiff_t)times.size(); ++t)
            {
                f << t << sep << times[t] << sep << states[t].cur_state.col(1ll).format(commaFmt) << '\n';
            }
            f.close();
        }

        // print grids
        {
            ofstream f{std::string{"output/z_grid.csv"}};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f.close();
        }
        {
            ofstream f{std::string{"output/r_grid.csv"}};
            f << grid_r.mesh_nodes.transpose().format(commaFmt) << '\n';
            f.close();
        }
        {
            ofstream f{std::string{"output/geotherma.csv"}};
            f << grid_z.mesh_nodes.transpose().format(commaFmt) << '\n';
            f << initial_state.cur_state.col(0ll).transpose().format(commaFmt) << '\n';
            f.close();
        }

        {
            ofstream f{std::string{"output/data.txt"}};
            f << "top collector height: " << grid_z.mesh_nodes(well.top_collector_cell_id()) << " m" << endl;
            f.close();
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
}

std::vector<RealType> Wrapper::get_times() const
{
    return time;
}
// std::vector<std::vector<RealType>> Wrapper::get_temps() const
// {
//     return t_radial_distribution;
// }

Wrapper::~Wrapper()
{
}