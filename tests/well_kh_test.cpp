#include <iostream>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Model/Well.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace std;

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::Phases;
using namespace GPN::CoordinateTypes;

using VR = std::vector<RealType>;

vector<RealType> make_steps(const VR &data)
{
    VR out(data.size() - 1ull);
    for (auto i{1ull}; i < data.size(); ++i)
        out[i - 1ull] = data[i] - data[i - 1ull];

    return out;
}

LogValuesContainer transfer_to_eigen(const VR &data, const RealType factor = 1.0)
{
    LogValuesContainer out(data.size());
    for (auto i{0ull}; i < data.size(); ++i)
        out(i) = factor * data[i];
    return out;
}

std::vector<RealType> grid_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
std::vector<RealType> grid_thickness{make_steps(grid_stencils)};
std::vector<RealType> permeability_stencils(grid_stencils.size() - 1ull, 1.0);
std::vector<RealType> is_permeable_stencils(grid_stencils.size() - 1ull, 1.0);
std::vector<RealType> is_perforated_stencils{is_permeable_stencils};

const RealType rate{1.0};

TEST_CASE("Well_KH_Test")
{
    is_perforated_stencils[0ll] = 0.0;

    const auto z_grid{
        Grids::Factory::create_axes<CoordinateTypes::Z>(
            grid_stencils)};

    const auto is_permeable{
        Logs::IsPermeableFactory::create(
            is_permeable_stencils,
            z_grid)};

    const auto is_perforated{
        Logs::IsPerforatedFactory::create(
            is_perforated_stencils,
            is_permeable_stencils,
            z_grid)};

    const auto permeability{
        Logs::PermeabilityFactory::create(
            permeability_stencils,
            is_permeable_stencils,
            z_grid)};

    const auto water{FluidFactory::create_water(1.0, 1.0)};

    const Well_KH well{
        water, is_permeable, is_perforated, permeability};

    const auto rfp = RFPFactory::create(well.get_RFP(rate), is_permeable);
    const auto wfp = WFPFactory::create(well.get_WFP(rate), is_perforated);

    cout << "thickness profile:     \n"
         << transfer_to_eigen(grid_thickness).transpose();

    cout << "\nwell rate profile:     \n"
         << wfp.log_vals.transpose();
    cout << "\nreservoir rate profile:\n"
         << rfp.log_vals.transpose();

    for (auto i{0ll}; i < rfp.size(); ++i)
    {
        CHECK(
            rfp(i) ==
            rate / (grid_stencils.back() - grid_stencils.front()) *
                z_grid.dual_steps(i));
    }

    CHECK(rfp.log_vals.sum() == rate);
    CHECK(wfp.log_vals.sum() == rate);

    RealType cum_rate{0.0};
    ptrdiff_t i{0ll};
    for (auto l{0ll}; (i < rfp.size()) && (l < 2ll); ++i)
    {
        if (is_permeable(i) == 1.0)
        {
            cum_rate += rfp(i);
            ++l;
        }
    }

    CHECK(wfp(i - 1ll) == cum_rate);
}
