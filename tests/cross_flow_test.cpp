#include <fstream>
#include <iostream>

#include <Injector/Grids/Defines.h>

#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Grids/GridRefiners.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/LogsFactory.hpp>

#include <Injector/Model/Well/CrossFlow.hpp>
#include <Injector/Model/Well/Well.hpp>

#include "includes/transfer_to_eigen.hpp"
#include "includes/set_is_permeable_stencils.hpp"

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;
using VR = std::vector<GPN::RealType>;

using namespace std;
using namespace Catch;
using namespace Catch::Matchers;
using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::CrossFlow;

RealType viscosity{6e-4}, density{1000}, capacity{4200}, heat_conductivity{0.6};

struct Well_CrossFlow
{
    Well_CrossFlow(
        const Logs::RFP &RFP_weights,
        const Logs::WFP &WFP_weights,
        const std::vector<RealType> &from_coords,
        const std::vector<ptrdiff_t> &to_layers)
        : RFP_weights{RFP_weights.log_vals},
          weights_sum{RFP_weights.log_vals.sum()},
          WFP_weights{WFP_weights.log_vals},
          cross_flows{RFP_weights, from_coords, to_layers}
    {
        assert(RFP_weights.log_vals.sum() == WFP_weights.log_vals.sum());
        assert((weights_sum == 1.0));
    }

    StepPropertyContainer get_RFP(
        RealType rate,
        RealType pressure) const
    {
        if (std::isnan(rate))
        { // define rate from pressure
            throw std::invalid_argument("RFP: Rate must be set");
        }
        else if (std::isnan(pressure))
        { // define pressure from rate
            assert(!std::isnan(rate));
            assert(rate >= 0.0);
            return ((rate / weights_sum) * RFP_weights).eval();
        }
        else
            throw std::invalid_argument("RFP: Either rate or pressure must be set, but not both.");
    }
    
    StepPropertyContainer get_WFP(
        RealType rate,
        RealType pressure) const
    {
        if (std::isnan(rate))
        { // define rate from pressure
            throw std::invalid_argument("RFP: Rate must be set");
        }
        else if (std::isnan(pressure))
        { // define pressure from rate
            assert(!std::isnan(rate));
            assert(rate >= 0.0);
            return ((rate / weights_sum) * WFP_weights).eval();
        }
        else
            throw std::invalid_argument("RFP: Either rate or pressure must be set, but not both.");
    }

private:
    const RealType weights_sum;
    const StepPropertyContainer &RFP_weights;
    const StepPropertyContainer &WFP_weights;

    const CrossFlows cross_flows;
};

TEST_CASE("CrossFlow", "")
{
    ifstream f("cross_flow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);

    const auto weights_stencils{transfer_to_eigen(data["collector"]["explicit"]["weights"].get<VR>())};

    const auto from_coords{data["collector"]["cross_flow"]["from_coord"].get<VR>()};
    const auto to_layers{data["collector"]["cross_flow"]["to_layers"].get<std::vector<std::ptrdiff_t>>()};
    const auto is_perforated_stencils{transfer_to_eigen(data["collector"]["is_perforated"].get<VR>())};
    const auto is_permeable_stencils{set_is_permeable_stencils(is_perforated_stencils, to_layers)};

    const auto grid_z{
        Factory::create_axes<CoordinateTypes::Z>(
            RefinerVerticle{
                data["grid"]["z_minor_step"].get<RealType>(),
                transfer_to_eigen(data["collector"]["is_permeable"].get<VR>())},
            Grids::Factory::generate_dual_grid_stencils_from_steps(
                0.0, data["collector"]["thickness"].get<VR>()))};

    const auto is_permeable{
        IsPermeableFactory::create(is_permeable_stencils, grid_z)};

    const auto weights{
        RFPFactory::create_from_container(
            weights_stencils,
            is_permeable)};

    const RealType total_rate{1.0};
    const CrossFlows cross_flows{
        weights, from_coords, to_layers};

    cout << "flux weights: " << weights.log_vals.transpose() << endl;
    const auto &dual_nodes{weights.grid.dual_nodes};
    const auto &dual_stencils{weights.grid.dual_stencils};
    for (auto i{0ull}; i < cross_flows.cross_flow_data.size(); ++i)
    {
        const auto &cf{cross_flows.cross_flow_data[i]};
        cout << "flux vector: " << cf.verticle_flux.transpose() << endl;

        auto to_id{0ll};
        { // get TO id
            const auto to_coord{dual_stencils(to_layers[i] + 1ll)};
            while (dual_nodes(to_id) < to_coord)
                ++to_id;
            --to_id;
            CHECK(dual_nodes(to_id) < to_coord);
            CHECK(dual_nodes(to_id + 1ll) == to_coord);
        }
        auto from_id{0ll};
        { // get FROM id
            while (dual_nodes(from_id) < from_coords[i])
                ++from_id;
            --from_id;
            CHECK(dual_nodes(from_id) < from_coords[i]);
            CHECK(dual_nodes(from_id + 1ll) > from_coords[i]);
        }

        const RealType dir{from_id < to_id ? 1.0 : -1.0};
        cout << "dir: " << dir << std::endl;
        for (auto j{0ll}; j < std::min(to_id, from_id) + 1ll; ++j)
        {
            INFO("top_id: " << j << ", to_id: " << to_id << ", from_id: " << from_id);
            CHECK(cf.verticle_flux(j) == 0.0);
        }
        for (auto j{std::max(to_id, from_id) + 1ll}; j < dual_nodes.rows(); ++j)
            CHECK(cf.verticle_flux(j) == 0.0);

        for (auto j{std::min(to_id, from_id) + 1ll}; j < std::max(to_id, from_id) + 1ll; ++j)
        {
            INFO("top_id: " << j << ", to_id: " << to_id << ", from_id: " << from_id << ", dir: " << dir);
            CHECK(cf.verticle_flux(j) == dir * weights(to_id));
        }
    }
}
