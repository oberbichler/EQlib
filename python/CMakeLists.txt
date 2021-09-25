cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(eqlib_python LANGUAGES CXX)

# --- Dependencies

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake")

find_package(MKL)

include(../cmake/CPM.cmake)

CPMAddPackage("gh:/pybind/pybind11@2.7.1")

CPMAddPackage("gh:gabime/spdlog@1.9.2")

CPMAddPackage("gh:Tessil/robin-map@0.6.3")

CPMAddPackage("gh:pybind/pybind11@2.7.1")

CPMAddPackage(
  NAME Eigen
  VERSION 3.4.0
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
  DOWNLOAD_ONLY YES 
)

if(Eigen_ADDED)
  add_library(Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()

if(TEST_INSTALLED_VERSION)
  find_package(eqlib REQUIRED)
else()
  CPMAddPackage(NAME eqlib SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# --- Create binary

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
pybind11_add_module(eqlib_python ${sources})

target_link_libraries(eqlib_python PRIVATE eqlib::eqlib spdlog Eigen tsl::robin_map MKL::RT)

set_target_properties(eqlib_python PROPERTIES CXX_STANDARD 17)

set_target_properties(eqlib_python PROPERTIES OUTPUT_NAME eqlib)