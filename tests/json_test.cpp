#include <fstream>
#include <nlohmann/json.hpp>
#include <catch2/catch_test_macros.hpp>


using json = nlohmann::json;
using namespace std;
// Tests Cylinder grid, (r; z)
TEST_CASE("BaseSplitTest")
{
    ifstream f("../../../tests/test_data/json_test_data.json");

    REQUIRE(f.is_open());
    json data = json::parse(f);

    CHECK(data["fluid"]["density"] == 1.0);
    CHECK(data["fluid"]["viscosity"] == 1.0);
    CHECK(data["fluid"]["specificHeatCapacity"] == 1.0);

    CHECK(data["rocks"]["nLayers"] == 10);
    CHECK(data["rocks"]["thickness"] == 1.0);
    CHECK(data["rocks"]["porosity"] == 1.0);
    CHECK(data["rocks"]["permeability"] == 1.0);
    CHECK(data["rocks"]["heatConductivity"] == 1.0);
    CHECK(data["rocks"]["solidDensity"] == 1.0);
    CHECK(data["rocks"]["soidSpecificHeatCapacity"] == 1.0);
    CHECK(data["rocks"]["initTemperature"] == 0.0);

    CHECK(data["grid"]["ztop"] == 0.0);
    CHECK(data["grid"]["z_step"] == 1.0);
    CHECK(data["grid"]["r_start"] == 1.0);
    CHECK(data["grid"]["r_end"] == 5.0);
    CHECK(data["grid"]["rnodes"] == 15);
    
    CHECK(data["history"]["t_start"] == 0.0);
    CHECK(data["history"]["t_end"] == 0.65);
    CHECK(data["history"]["t_step"] == 0.1);
    CHECK(data["history"]["wellRate"] == 2.0);
    CHECK(data["history"]["inletTemp"] == 1.0);
}
