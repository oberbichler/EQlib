cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)

project(EQlib-py)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(TBB REQUIRED)

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
    message("-- Found MKL: $ENV{CONDA_PREFIX}/Library")
    message("-- Found TBB: $ENV{CONDA_PREFIX}/Library")

    include_directories(INTERFACE $ENV{CONDA_PREFIX}/Library/include)
    link_directories($ENV{CONDA_PREFIX}/Library/lib)

    add_definitions(-DEIGEN_USE_MKL_ALL)

    target_link_libraries(EQlib PRIVATE mkl_rt)
    target_link_libraries(EQlib PRIVATE tbb.lib debug tbb_debug.lib)
else()
    if(EXISTS $ENV{MKLROOT})
        message("-- Found MKL: $ENV{MKLROOT}")
        include_directories(INTERFACE $ENV{MKLROOT}/include)
        add_definitions(-DEIGEN_USE_MKL_ALL)
        target_link_libraries(EQlib PRIVATE mkl_rt)
    endif()

    if(${TBB_FOUND})
        message("-- Found TBB: $ENV{TBBROOT}")
        include_directories(INTERFACE ${TBB_INCLUDE_DIRS})
        target_compile_definitions(EQlib PRIVATE ${TBB_DEFINITIONS})
        target_link_libraries(EQlib PRIVATE ${TBB_LIBRARIES})
    endif()
endif()

install(TARGETS EQlib DESTINATION bin)
install(TARGETS IGAlib DESTINATION bin)