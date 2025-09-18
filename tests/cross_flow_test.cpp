#include <fstream>

#include <Injector/Grids/Defines.h>

#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;
using VR = std::vector<GPN::RealType>;

using namespace std;
using namespace Catch;
using namespace Catch::Matchers;
using namespace GPN;

RealType viscosity{6e-4}, density{1000}, capacity{4200}, heat_conductivity{0.6};

TEST_CASE("CrossFlow", "Water")
{
    ifstream f("cross_flow_test_data.json");
    REQUIRE(f.is_open());
    json data = json::parse(f);



}
