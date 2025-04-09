# https://google.github.io/googletest/primer.html
# https://habr.com/ru/articles/667880/
# https://habr.com/ru/articles/119090/

cmake_minimum_required(VERSION 3.14)

# create a variable called target_name and set it to the string "dtest"
set (target_name injector_test)
project(${target_name})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/..") 


include(CTest)

macro(add_google_test name)
  set(TEST_EXECUTABLE ${name})
  set(TEST_NAME ${name})

  add_executable(
    ${TEST_NAME}  
      ${TEST_EXECUTABLE}.cpp
  )
  target_include_directories(
    ${TEST_NAME} PRIVATE 
      "${CMAKE_CURRENT_SOURCE_DIR}/../"
  )
  target_include_directories(
    ${TEST_NAME} PRIVATE 
      "${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen"
  )
  target_include_directories(
    ${TEST_NAME} PRIVATE 
      "${CMAKE_CURRENT_SOURCE_DIR}/../deps/gtest"
  )
  target_link_libraries(
    ${TEST_NAME} PRIVATE
      gtest_main
  ) 

  # let the preprocessor know about the system name
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_compile_definitions(${TEST_NAME} PUBLIC "IS_LINUX")
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    target_compile_definitions(${TEST_NAME} PUBLIC "IS_MACOS")
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    target_compile_definitions(${TEST_NAME} PUBLIC "IS_WINDOWS")
  endif()

  # compiler-specific options
  if(MSVC AND MSVC_VERSION GREATER 1400)
    target_compile_options(
      ${TEST_NAME} PRIVATE 
      /MP)
  else()
    target_compile_options(
      ${TEST_NAME} PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
  endif()
  
  if(MSVC)
    # hack for windows compilation
    # or for compilation with MSVC
    # if MinGW -- may not be required...
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()

  # adding the test, finally
  add_test(
    NAME ${TEST_NAME} 
    COMMAND ${TEST_EXECUTABLE})

endmacro()

# add_google_test(point_test)
# add_google_test(pointval_test)


