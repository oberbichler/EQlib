cmake_minimum_required(VERSION 3.12)

# get commit hash
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(CMAKE_CXX_STANDARD 17)

project(eqlib-py)

include_directories(
    "include"
    ${CMAKE_BINARY_DIR}/generated
    SYSTEM "external_libraries"
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(MKL)

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
    -DSPDLOG_FMT_EXTERNAL
)

pybind11_add_module(eqlib src/Module.cpp)

if (DEFINED ENV{CONDA_PREFIX})
    message("-- Found Anaconda: $ENV{CONDA_PREFIX}")

    if(WIN32)
        include_directories(SYSTEM $ENV{CONDA_PREFIX}/Library/include)
        set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "$ENV{CONDA_PREFIX}/Library/lib")
    else()
        include_directories(SYSTEM $ENV{CONDA_PREFIX}/include)
        set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "$ENV{CONDA_PREFIX}/lib")
    endif()
else()
    message("-- Anaconda not found")
endif()

if (MKL_FOUND)
    message("-- Found MKL: ${MKL_INCLUDE_DIR}")
    include_directories(SYSTEM ${MKL_INCLUDE_DIR})
    target_link_libraries(eqlib PRIVATE ${MKL_LIBRARIES})
else()
    message("-- MKL not found")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(eqlib PUBLIC -fsized-deallocation)
endif()


configure_file(
  ${CMAKE_SOURCE_DIR}/include/eqlib/Info.h.in
  ${CMAKE_BINARY_DIR}/generated/eqlib/Info.h
)


install(TARGETS eqlib DESTINATION bin)
