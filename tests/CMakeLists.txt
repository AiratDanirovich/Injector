# https://google.github.io/googletest/primer.html
# https://habr.com/ru/articles/667880/
# https://habr.com/ru/articles/119090/

cmake_minimum_required(VERSION 3.14)

# create a variable called project_name and set it to the string "injector_test"
set (project_name injector_test)
project(${project_name})

# Include Catch2
include(${CMAKE_SOURCE_DIR}../cmake/Catch2.cmake)
# add macro with add_catch2_test definition
include(${CMAKE_SOURCE_DIR}../cmake/AddTest.cmake)

# Enable CTest
include(CTest)
include(Catch)
# catch_discover_tests(tests)

# include_directories("${CMAKE_CURRENT_SOURCE_DIR}/..") 

# Create test executable
add_catch2_test(water_test)
add_catch2_test(stride_test)
add_catch2_test(grid_test)
add_catch2_test(initialstate_test)
add_catch2_test(boundaryconditions_test)
add_catch2_test(logs_test)
add_catch2_test(physicalfield_test)
add_catch2_test(history_test)
add_catch2_test(basesplit_test)
add_catch2_test(splitx_test)
add_catch2_test(splity_test)
add_catch2_test(solver_test)
add_catch2_test(splitting_solver_constvalProblem)
add_catch2_test(splitting_solver_SelfSimilarProblem)


