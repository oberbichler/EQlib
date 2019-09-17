cmake_minimum_required(VERSION 3.12)

# get commit hash
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(CMAKE_CXX_STANDARD 20)

project(EQlib-py)

configure_file(
  ${CMAKE_SOURCE_DIR}/include/EQlib/Version.h.in
  ${CMAKE_BINARY_DIR}/generated/EQlib/Version.h
)

include_directories(
    "include"
    ${CMAKE_BINARY_DIR}/generated
    SYSTEM "external_libraries"
)

find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(external_libraries/pybind11)

add_definitions(
    -DEIGEN_DEFAULT_TO_ROW_MAJOR
    -DFMT_HEADER_ONLY
)

pybind11_add_module(EQlib src/module.cpp)
pybind11_add_module(IGAlib src/iga_module.cpp)

if (DEFINED ENV{CONDA_PREFIX})
    message("-- Found Anaconda: $ENV{CONDA_PREFIX}")

    if(WIN32)
        include_directories(SYSTEM $ENV{CONDA_PREFIX}/Library/include)
    else()
        include_directories(SYSTEM $ENV{CONDA_PREFIX}/include)
    endif()
else()
    message("-- Anaconda not found")
endif()


find_library(MKL_LIBRARY mkl_rt)

if(MKL_LIBRARY)
    message("-- Found MKL: ${MKL_LIBRARY}")
    add_definitions(-DEIGEN_USE_MKL_ALL)
    target_link_libraries(EQlib PRIVATE ${MKL_LIBRARY})
    target_link_libraries(IGAlib PRIVATE ${MKL_LIBRARY})
else()
    message(FATAL_ERROR "-- MKL not found")
endif()


find_library(TBB_LIBRARY tbb)
find_library(TBB_DBG_LIBRARY tbb_debug)

if(TBB_LIBRARY)
    message("-- Found TBB: ${TBB_LIBRARY}")
    target_link_libraries(EQlib PRIVATE ${TBB_LIBRARY} debug ${TBB_DBG_LIBRARY})
    target_link_libraries(IGAlib PRIVATE ${TBB_LIBRARY} debug ${TBB_DBG_LIBRARY})
else()
    message(FATAL_ERROR "-- TBB not found")
endif()


if(DEFINED WORHP_ROOT)
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${WORHP_ROOT}/bin")
endif()

find_library(WORHP_LIBRARY worhp)

if(WORHP_LIBRARY)
    message("-- Found WORHP: ${WORHP_LIBRARY}")
    add_definitions(-DUSE_WORHP)
    target_include_directories(EQlib SYSTEM PRIVATE "${WORHP_ROOT}/include")
    target_link_libraries(EQlib PRIVATE ${WORHP_LIBRARY})
else()
    message("-- WORHP not found")
endif()


install(TARGETS EQlib IGAlib DESTINATION bin)
