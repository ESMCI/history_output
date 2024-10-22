cmake_minimum_required(VERSION 3.21)

project(History-Output VERSION 0.0.1 LANGUAGES Fortran)

set(CMAKE_Fortran_FLAGS "-O0 --coverage")
if(NOT History-Output_IS_TOP_LEVEL)
  message(WARNING "History-Output is not integrated into the CMake build of any top level "
                  "project yet and this CMake is for testing purposes only.  "
                  "Making a change to this project's CMake will not impact the build of "
                  "a parent project at this time.")
endif()

option(HISTORY_ENABLE_TESTS "Run pFUnit unit tests" OFF)
option(HISTORY_ENABLE_CODE_COVERAGE "Run code coverage tool" OFF)

if(HISTORY_ENABLE_CODE_COVERAGE)
  add_compile_options(-O0 --coverage)
  add_link_options(-lgcov)
endif()

set(DCMAKE_BUILD_TYPE Debug)
add_subdirectory(src)

if(HISTORY_ENABLE_TESTS OR HISTORY_ENABLE_CODE_COVERAGE)
  enable_testing()
  add_subdirectory(test)
endif()
