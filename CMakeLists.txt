# useful links
# https://cgold.readthedocs.io/en/latest/tutorials/libraries/static-shared.html
# https://cgold.readthedocs.io/en/latest/overview/cmake-can.html
# https://stackoverflow.com/questions/17511496/how-to-create-a-shared-library-with-cmake
# https://vector-of-bool.github.io/docs/vscode-cmake-tools/configuring.html
# https://coderefinery.github.io/cmake-workshop/flags-definitions-debugging/
# https://cmake.org/cmake/help/latest/guide/tutorial/Packaging%20Debug%20and%20Release.html
# https://coderefinery.github.io/cmake-workshop/fetch-content/

# https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
# https://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake
# https://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake

# https://github.com/esufaev/Parallel-Operations-Toolkit
# https://github.com/Data-Analysis-Department/temp_calc/blob/py/build.cmd

# https://ps-group.github.io/cxx/cmake_cheatsheet
# https://batuhankoc.medium.com/use-cmake-to-create-and-use-dynamic-libraries-5f6498417b3c

# https://coderefinery.github.io/cmake-workshop/
# https://cliutils.gitlab.io/modern-cmake/
# https://cmake.org/cmake/help/v3.19/index.html

# to add submodule
# git submodule add https://github.com/google/googletest deps/gtest/
# git submodule add https://gitlab.com/libeigen/eigen deps/eigen/

cmake_minimum_required(VERSION 3.14)

project(InjectorGPN VERSION 1.0 LANGUAGES C CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

find_package(OpenMP)
    if (OPENMP_FOUND)
            message(STATUS "OpenMP Found") 
            set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
            list(APPEND GPN_USE_LIST OPENMP) # GPN_USE_LIST may be required...
    else()
            message(STATUS "OpenMP NOT Found") 
    endif()

# enable testing
enable_testing()
# process tests/CMakeLists.txt
add_subdirectory(tests)

add_subdirectory(InjectorDLL)
add_subdirectory(InjectorUser)


install(TARGETS InjectorUser
    # copy *.exe file
    DESTINATION ${PROJECT_BINARY_DIR}/bin
    )
install(TARGETS InjectorGPNDLL
    # copy *.dll file
    RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin
    # copy public header(s)
    PUBLIC_HEADER DESTINATION ${PROJECT_BINARY_DIR}/bin
    )
    # message(PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR})
install(FILES ${PROJECT_SOURCE_DIR}/tests/test_data/injector_launch.json
    DESTINATION ${PROJECT_BINARY_DIR}/bin)


# install section
# to copy generated binaries
# install(TARGETS dll_user
#     # copy *.exe file
#     DESTINATION ${PROJECT_BINARY_DIR}/bin
#     )
    
# install(TARGETS TestDLL
#     # copy *.dll file
#     RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin
#     # copy public header(s)
#     PUBLIC_HEADER DESTINATION ${PROJECT_BINARY_DIR}/bin
#     )
# install(FILES ${PROJECT_BINARY_DIR}/dll_src/Debug/TestDLL.lib
#     # copy *.lib file
#     DESTINATION ${PROJECT_BINARY_DIR}/bin)