cmake_minimum_required(VERSION 3.14)
project(InjectorGPNPyd)

set(CMAKE_CXX_STANDARD 17)

# Указать вручную путь к pybind11
include_directories("C:/Users/KFU_USER_1/AppData/Roaming/Python/Python311/site-packages/pybind11/include")
# Путь к include каталогу Python
include_directories("C:/ProgramData/anaconda3/include")
# Путь к .lib файлам Python
link_directories("C:/ProgramData/anaconda3/libs")

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_SOURCE_DIR})

add_library(InjectorGPN MODULE bindings.cpp)

set_target_properties(InjectorGPN PROPERTIES PREFIX "")
set_target_properties(InjectorGPN PROPERTIES SUFFIX ".pyd")

target_link_libraries(InjectorGPN PRIVATE InjectorGPNDLL)
