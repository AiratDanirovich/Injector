add_executable(InjectorUser main.cpp)
target_link_libraries(InjectorUser PUBLIC InjectorGPNDLL)
target_link_libraries(InjectorUser PUBLIC nlohmann_json::nlohmann_json)

target_include_directories(InjectorUser PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          )