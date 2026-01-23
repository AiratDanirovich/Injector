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

    cout << "Simulation is started." << endl;
    cout << "Please wait..." << endl;
    const auto t_start{chrono::high_resolution_clock::now()};
    Wrapper *instance = new Wrapper(data);
    const auto t_end{chrono::high_resolution_clock::now()};
    cout << "Elapsed time:                       " << (t_end - t_start).count() * 1E-9 << " seconds\n";

    const auto &t = instance->get_times();
    cout << "number of saved time moments:       " << t.size() << endl;

    delete instance;

    std::cout << "Simulation is finished.\nPress any key to exit..." << std::endl;
    getchar();

    return 0;
}
