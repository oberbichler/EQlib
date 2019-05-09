cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)

project(EQlib-py)

include_directories(
    "include"
    SYSTEM "external_libraries"
)

add_subdirectory(external_libraries/pybind11)

add_definitions(
    -DEIGEN_DEFAULT_TO_ROW_MAJOR
    -DEQlib_VERSION="${EQlib_VERSION}"
    -DFMT_HEADER_ONLY
)

pybind11_add_module(EQlib src/module.cpp)
pybind11_add_module(IGAlib src/iga_module.cpp)

if (DEFINED ENV{CONDA_PREFIX})
    message("-- Found Anaconda: $ENV{CONDA_PREFIX}")

    include_directories(SYSTEM $ENV{CONDA_PREFIX}/Library/include)
    include_directories(SYSTEM $ENV{CONDA_PREFIX}/include)

    add_definitions(-DEIGEN_USE_MKL_ALL)

    find_library(MKL_LIBRARY mkl_rt HINTS $ENV{CONDA_PREFIX}/Library/lib $ENV{CONDA_PREFIX}/lib)
    find_library(TBB_LIBRARY tbb HINTS $ENV{CONDA_PREFIX}/Library/lib $ENV{CONDA_PREFIX}/lib)
    find_library(TBB_DBG_LIBRARY tbb_debug HINTS $ENV{CONDA_PREFIX}/Library/lib $ENV{CONDA_PREFIX}/lib)

    target_link_libraries(EQlib PRIVATE ${MKL_LIBRARY})
    target_link_libraries(EQlib PRIVATE ${TBB_LIBRARY} debug ${TBB_DBG_LIBRARY})

    target_link_libraries(IGAlib PRIVATE ${MKL_LIBRARY})
    target_link_libraries(IGAlib PRIVATE ${TBB_LIBRARY} debug ${TBB_DBG_LIBRARY})
endif()

install(TARGETS EQlib DESTINATION bin)
install(TARGETS IGAlib DESTINATION bin)
