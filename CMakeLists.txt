cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)

project(EQLib-py)

find_package(BLAS)

include_directories(
    "include"
    "external_libraries"
)

add_subdirectory(external_libraries/pybind11)

add_definitions(
    -DEIGEN_DEFAULT_TO_ROW_MAJOR
    -DEQLIB_VERSION="${EQLIB_VERSION}"
    -DFMT_HEADER_ONLY
)

pybind11_add_module(EQLib src/module.cpp)

if (EXISTS $ENV{MKLROOT})
    message("Found MKL: $ENV{MKLROOT}")
    include_directories($ENV{MKLROOT}/include)
    add_definitions(-DEIGEN_USE_MKL_ALL)
    target_link_libraries(EQLib PRIVATE mkl_rt)
elseif(BLAS_FOUND)
    add_definitions(-DEIGEN_USE_BLAS)
    target_link_libraries(EQLib PRIVATE ${BLAS_LIBRARIES})
endif()

install(TARGETS EQLib DESTINATION bin)