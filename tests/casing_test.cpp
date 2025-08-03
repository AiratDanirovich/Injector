#include <fstream>
#include <limits>

#include <Injector/Model/Completion.hpp>
#include <Injector/Model/CasingLogs.hpp>

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace std;

using namespace Catch;
using namespace Catch::Matchers;
using json = nlohmann::json;

using namespace GPN;
using namespace GPN::Completion;

TEST_CASE("Well_Test")
{
    ifstream f("heatflow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);

    vector<VarRing> casing;

    { // flowing fluid
        const auto &data2 = data["fluid"];
        casing.push_back(
            VarRing{
                VarFlow{
                    Completion::Density{data2["density"].get<vector<RealType>>()},
                    Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                    Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
                Completion::VarThickness{data["completion"]["tube"]["inner_radius"].get<vector<RealType>>()},
                Completion::VarInnerRadius{vector{0.0}},
                Completion::VarDepth{vector{std::numeric_limits<RealType>::max()}}});
    }
    {
    const auto &data2 = data["completion"]["variable"]["tube"];
    casing.emplace_back(
        VarRing {
            VarTube{
                Completion::Density{data2["density"].get<vector<RealType>>()},
                Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
            Completion::VarThickness{data2["thickness"].get<vector<RealType>>()},
            Completion::VarInnerRadius{casing.back().outer_radius},
            Completion::VarDepth{data2["heat_conductivity"].get<vector<RealType>>()} });
    }    
    {
    const auto &data2 = data["completion"]["variable"]["annulus"];
    casing.emplace_back(
        VarRing {
            VarAnnulus{
                Completion::Density{data2["density"].get<vector<RealType>>()},
                Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
            Completion::VarThickness{data2["thickness"].get<vector<RealType>>()},
            Completion::VarInnerRadius{casing.back().outer_radius},
            Completion::VarDepth{data2["heat_conductivity"].get<vector<RealType>>()} });
    }    
    {
    const auto &data2 = data["completion"]["variable"]["column"];
    casing.emplace_back(
        VarRing {
            VarColumn{
                Completion::Density{data2["density"].get<vector<RealType>>()},
                Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
            Completion::VarThickness{data2["thickness"].get<vector<RealType>>()},
            Completion::VarInnerRadius{casing.back().outer_radius},
            Completion::VarDepth{data2["heat_conductivity"].get<vector<RealType>>()} });
    }    {
    const auto &data2 = data["completion"]["variable"]["cement"];
    casing.emplace_back(
        VarRing {
            VarCement{
                Completion::Density{data2["density"].get<vector<RealType>>()},
                Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
            Completion::VarThickness{data2["thickness"].get<vector<RealType>>()},
            Completion::VarInnerRadius{casing.back().outer_radius},
            Completion::VarDepth{data2["heat_conductivity"].get<vector<RealType>>()} });
    }


    Casing<VarRing> completion{casing};
}
