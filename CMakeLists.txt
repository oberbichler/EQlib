cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# --- Project

project(
    eqlib
    VERSION 2.0.0
    LANGUAGES CXX
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- Include guards

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
    )
endif()

# --- Custom CMake

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# --- Get CPM

include(cmake/CPM.cmake)

# --- Add dependencies

find_package(MKL)

include(cmake/CPM.cmake)

CPMAddPackage(
    NAME HyperJet
    VERSION 1.4.3
    URL https://github.com/oberbichler/HyperJet/archive/refs/tags/v1.4.3.tar.gz
    DOWNLOAD_ONLY YES 
)

CPMAddPackage("gh:gabime/spdlog@1.10.0")
set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)

CPMAddPackage("gh:Tessil/robin-map@1.0.1")

CPMAddPackage("gh:pybind/pybind11@2.9.2")

CPMAddPackage(
    NAME Eigen
    VERSION 3.4.0
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_ONLY YES 
)

if(Eigen_ADDED)
    add_library(eigen::eigen INTERFACE IMPORTED)
    target_include_directories(eigen::eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()

# ---

add_subdirectory(${pybind11_SOURCE_DIR})

pybind11_add_module(python_module
    # eqlib::
    src/python/main.cpp
    src/python/constraint.cpp
    src/python/problem.cpp
    src/python/equation.cpp
    src/python/linear_solver.cpp
    src/python/node.cpp
    src/python/objective.cpp
    src/python/parameter.cpp
    src/python/request.cpp
    src/python/sparse_structure.cpp
    # eqlib::objectives
    src/python/objectives/lambda_objective.cpp
    src/python/objectives/spring.cpp
    # eqlib::constraints
    src/python/constraints/lambda_constraint.cpp
    # eqlib::linear_solvers
    src/python/linear_solvers/pardiso_ldlt.cpp
    src/python/linear_solvers/simplicial_ldlt.cpp
)

target_compile_features(python_module PRIVATE cxx_std_17)

target_include_directories(python_module PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(MSVC)
    set(EQLIB_OUTPUT_DIRECTORY
        LIBRARY_OUTPUT_DIRECTORY_RELEASE        ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
        LIBRARY_OUTPUT_DIRECTORY_DEBUG          ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
        LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
        LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL     ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
    )
else()
    set(EQLIB_OUTPUT_DIRECTORY
        LIBRARY_OUTPUT_DIRECTORY                ${CMAKE_CURRENT_BINARY_DIR}/../../eqlib
    )
endif()

target_link_libraries(python_module PRIVATE
    eigen::eigen
    spdlog::spdlog
    tsl::robin_map
)

if(${MKL_FOUND})
    target_compile_definitions(python_module PRIVATE EIGEN_USE_MKL_ALL EQLIB_USE_MKL)
    target_link_libraries(python_module PRIVATE MKL::MKL)
endif()

set_target_properties(python_module PROPERTIES OUTPUT_NAME eqlib_ext ${EQLIB_OUTPUT_DIRECTORY})

install(TARGETS python_module LIBRARY DESTINATION eqlib)