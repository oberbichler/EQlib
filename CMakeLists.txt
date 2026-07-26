cmake_minimum_required(VERSION 3.24)

project(eqlib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- version info

if(NOT DEFINED EQLIB_VERSION)
    if(DEFINED SKBUILD_PROJECT_VERSION)
        set(EQLIB_VERSION ${SKBUILD_PROJECT_VERSION})
    else()
        set(EQLIB_VERSION "dev")
    endif()
endif()

execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/eqlib/Info.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/eqlib/Info.h
)

# --- dependencies

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# header-only dependencies: only the sources are downloaded, their build
# systems are skipped (SOURCE_SUBDIR points to a directory without CMakeLists)
FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    SOURCE_SUBDIR cmake_disabled
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_Declare(fmt
    URL https://github.com/fmtlib/fmt/archive/refs/tags/11.1.4.tar.gz
    SOURCE_SUBDIR cmake_disabled
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_Declare(spdlog
    URL https://github.com/gabime/spdlog/archive/refs/tags/v1.15.2.tar.gz
    SOURCE_SUBDIR cmake_disabled
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_Declare(robin_map
    URL https://github.com/Tessil/robin-map/archive/refs/tags/v1.3.0.tar.gz
    SOURCE_SUBDIR cmake_disabled
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(eigen fmt spdlog robin_map)

# --- module

pybind11_add_module(eqlib
    src/Module.cpp
    src/bindings/core.cpp
    src/bindings/elements.cpp
    src/bindings/objectives.cpp
    src/bindings/problem.cpp
    src/bindings/solvers.cpp
)

target_include_directories(eqlib PRIVATE
    include
    ${CMAKE_CURRENT_BINARY_DIR}/generated
)

target_include_directories(eqlib SYSTEM PRIVATE
    ${eigen_SOURCE_DIR}
    ${fmt_SOURCE_DIR}/include
    ${spdlog_SOURCE_DIR}/include
    ${robin_map_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries
)

target_compile_definitions(eqlib PRIVATE
    EIGEN_DEFAULT_TO_ROW_MAJOR
    FMT_HEADER_ONLY
    SPDLOG_FMT_EXTERNAL
)

if(MSVC)
    # fmt >= 11 requires the UTF-8 source and execution charset
    target_compile_options(eqlib PRIVATE /utf-8)
endif()

# --- OpenMP (optional)

find_package(OpenMP COMPONENTS CXX)

if(OpenMP_CXX_FOUND)
    target_link_libraries(eqlib PRIVATE OpenMP::OpenMP_CXX)
endif()

# --- MKL (optional)

if(DEFINED ENV{MKL_ROOT})
    target_include_directories(eqlib SYSTEM PRIVATE $ENV{MKL_ROOT}/include)
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "$ENV{MKL_ROOT}/lib")
endif()

find_library(MKL_LIBRARY mkl_rt)

if(MKL_LIBRARY)
    message(STATUS "Found MKL: ${MKL_LIBRARY}")
    target_compile_definitions(eqlib PRIVATE EIGEN_USE_MKL_ALL EQLIB_USE_MKL)
    target_link_libraries(eqlib PRIVATE ${MKL_LIBRARY})
endif()

# --- install

install(TARGETS eqlib LIBRARY DESTINATION .)
