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

    const auto &data2 = data["completion"]["var_cement"];

    const auto cement_log{
        VarRing{
            VarCement{
                Completion::Density{data2["density"].get<vector<RealType>>()},
                Completion::SpecificHeatCapacity{data2["specific_heat_capacity"].get<vector<RealType>>()},
                Completion::HeatConductivity{data2["heat_conductivity"].get<vector<RealType>>()}},
            Completion::VarThickness{data2["thickness"].get<vector<RealType>>()},
            Completion::VarInnerRadius{data2["specific_heat_capacity"].get<vector<RealType>>()},
            Completion::VarDepth{data2["heat_conductivity"].get<vector<RealType>>()}}};


    
}
