#include <fstream>
#include <limits>

#include <Injector/Model/Completion.hpp>

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

    const auto &data2 = data["completion"]["cement"];

    const auto cement_log{
        Ring{
            Cement{
                Density{data2["density"]},
                SpecificHeatCapacity{data2["specific_heat_capacity"]},
                GPN::HeatConductivity{data2["heat_conductivity"]}},
            Thickness{data2["thickness"]},
            InnerRadius{5.0},
            Depth{std::numeric_limits<RealType>::max()}}
    };
}
