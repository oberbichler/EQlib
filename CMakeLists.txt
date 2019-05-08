cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)

project(EQlib-py)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include_directories(
    "include"
    INTERFACE "external_libraries"
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
    message("-- Found MKL: $ENV{CONDA_PREFIX}")
    message("-- Found TBB: $ENV{CONDA_PREFIX}")

    include_directories(INTERFACE $ENV{CONDA_PREFIX}/include)

    add_definitions(-DEIGEN_USE_MKL_ALL)

    find_library(MKL_LIBRARY mkl_rt HINTS $ENV{CONDA_PREFIX}/lib)
    find_library(TBB_LIBRARY tbb HINTS $ENV{CONDA_PREFIX}/lib)
    find_library(TBB_DBG_LIBRARY tbb_debug HINTS $ENV{CONDA_PREFIX}/lib)

    target_link_libraries(EQlib PRIVATE ${MKL_LIBRARY})
    target_link_libraries(EQlib PRIVATE ${TBB_LIBRARY} debug ${TBB_DBG_LIBRARY})
else()
    if(EXISTS $ENV{MKLROOT})
        message("-- Found MKL: $ENV{MKLROOT}")
        include_directories(INTERFACE $ENV{MKLROOT}/include)
        add_definitions(-DEIGEN_USE_MKL_ALL)
        target_link_libraries(EQlib PRIVATE mkl_rt)
    endif()

    find_package(TBB REQUIRED)

    if(${TBB_FOUND})
        message("-- Found TBB: $ENV{TBBROOT}")
        include_directories(INTERFACE ${TBB_INCLUDE_DIRS})
        target_compile_definitions(EQlib PRIVATE ${TBB_DEFINITIONS})
        target_link_libraries(EQlib PRIVATE ${TBB_LIBRARIES})
    endif()
endif()

install(TARGETS EQlib DESTINATION bin)
install(TARGETS IGAlib DESTINATION bin)