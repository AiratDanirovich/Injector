# project(InjectorGPNDLL VERSION 1.0 LANGUAGES C CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# find_package(OpenMP)
# if(OPENMP_FOUND)
#     message(STATUS "OpenMP Found")
#     if(NOT ${TBFMM_ENABLE_OPENMP})
#         message(STATUS "OpenMP Disabled")
#     else()
#         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#         set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#         list(APPEND GPN_USE_LIST OPENMP) # GPN_USE_LIST may be required...
#     endif()
# endif()

# Add the source files to the project
set(SOURCE_FILES
    Wrapper.cpp
)

# Add the header files to the project
set(HEADER_FILES
    Wrapper.h
)

if (MSVC)
    # error C1128: number of sections exceeded object file format limit: compile with /bigobj
    add_compile_options(/bigobj)
endif()

# Create the shared library (DLL)
add_library(InjectorGPNDLL SHARED ${HEADER_FILES}
    ${SOURCE_FILES})
target_link_libraries(InjectorGPNDLL PUBLIC nlohmann_json::nlohmann_json)
set_target_properties(InjectorGPNDLL PROPERTIES
    VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER InjectorDLL/Wrapper.h
)
target_compile_definitions(InjectorGPNDLL PRIVATE MYLIBRARY_EXPORT)
target_include_directories(InjectorGPNDLL INTERFACE
    "${PROJECT_SOURCE_DIR}"
)
target_include_directories(InjectorGPNDLL PRIVATE
    "${PROJECT_SOURCE_DIR}"
)
# directory with Eigen headers
target_include_directories(
    InjectorGPNDLL PRIVATE 
      "${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen"
  )

message(PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR})


