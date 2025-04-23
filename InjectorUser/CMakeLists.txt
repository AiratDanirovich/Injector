add_executable(InjectorUser main.cpp)
target_link_libraries(InjectorUser PUBLIC InjectorGPNDLL)

target_include_directories(InjectorUser PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          )