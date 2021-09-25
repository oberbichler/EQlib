cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# ---- Project ----

project(
  eqlib
  VERSION 2.0.0
  LANGUAGES CXX
)

#

set(CPM_DOWNLOAD_VERSION 0.32.0)

if(CPM_SOURCE_CACHE)
  set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
elseif(DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
else()
  set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
endif()

if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
  message(STATUS "Downloading CPM.cmake to ${CPM_DOWNLOAD_LOCATION}")
  file(DOWNLOAD
       https://github.com/TheLartians/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
       ${CPM_DOWNLOAD_LOCATION}
  )
endif()

include(${CPM_DOWNLOAD_LOCATION})

# ---- Add dependencies ----

include(cmake/CPM.cmake)

CPMAddPackage(
  NAME HyperJet
  VERSION 1.4.3
  URL https://github.com/oberbichler/HyperJet/archive/refs/tags/v1.4.3.tar.gz
  DOWNLOAD_ONLY YES 
)

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

# ---

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(MKL REQUIRED)

# ---

add_subdirectory(${pybind11_SOURCE_DIR})

pybind11_add_module(python_module
  src/python/main.cpp
  src/python/constraint.cpp
  src/python/problem.cpp
  src/python/equation.cpp
  src/python/node.cpp
  src/python/objective.cpp
  src/python/request.cpp
  src/python/sparse_structure.cpp
  src/python/variable.cpp
  # objectives
  src/python/lambda_constraint.cpp
  src/python/lambda_objective.cpp
  src/python/objectives/spring.cpp
)

target_compile_features(python_module PRIVATE cxx_std_17)

target_compile_definitions(python_module PRIVATE
  -DEQLIB_VERSION="${EQLIB_VERSION}"
)

target_include_directories(python_module PRIVATE
  include
  external/robin-map/include
  external/spdlog/include
)

if(MSVC)
  set(EQLIB_OUTPUT_DIRECTORY
    LIBRARY_OUTPUT_DIRECTORY_RELEASE        ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
    LIBRARY_OUTPUT_DIRECTORY_DEBUG          ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
    LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
    LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL     ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib)
else()
  set(EQLIB_OUTPUT_DIRECTORY
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib)
endif()

target_link_libraries(python_module PRIVATE spdlog Eigen tsl::robin_map MKL::RT)

set_target_properties(python_module PROPERTIES OUTPUT_NAME eqlib_ext ${EQLIB_OUTPUT_DIRECTORY})

install(TARGETS python_module LIBRARY DESTINATION eqlib)

# install(TARGETS python_module DESTINATION bin)