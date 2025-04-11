# useful links
# https://cgold.readthedocs.io/en/latest/tutorials/libraries/static-shared.html
# https://cgold.readthedocs.io/en/latest/overview/cmake-can.html
# https://stackoverflow.com/questions/17511496/how-to-create-a-shared-library-with-cmake

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

# enable testing
enable_testing()
# process tests/CMakeLists.txt
add_subdirectory(tests)


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