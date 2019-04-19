cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)

project(EQlib-py)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(BLAS)
find_package(TBB REQUIRED)

include_directories(
    "include"
    "external_libraries"
    "${HYPERJET_ROOT}/include"
)

add_subdirectory(external_libraries/pybind11)

add_definitions(
    -DEIGEN_DEFAULT_TO_ROW_MAJOR
    -DEQlib_VERSION="${EQlib_VERSION}"
    -DFMT_HEADER_ONLY
)

pybind11_add_module(EQlib src/module.cpp)
pybind11_add_module(IGAlib src/iga_module.cpp)

if(EXISTS $ENV{MKLROOT})
    message("-- Found MKL: $ENV{MKLROOT}")
    include_directories($ENV{MKLROOT}/include)
    add_definitions(-DEIGEN_USE_MKL_ALL)
    target_link_libraries(EQlib PRIVATE mkl_rt)
elseif(BLAS_FOUND)
    add_definitions(-DEIGEN_USE_BLAS)
    target_link_libraries(EQlib PRIVATE ${BLAS_LIBRARIES})
endif()

if(${TBB_FOUND})
    message("-- Found TBB: $ENV{TBBROOT}")
    include_directories(${TBB_INCLUDE_DIRS})
    target_compile_definitions(EQlib PRIVATE ${TBB_DEFINITIONS})
    target_link_libraries(EQlib PRIVATE ${TBB_LIBRARIES})
endif()

install(TARGETS EQlib DESTINATION bin)
install(TARGETS IGAlib DESTINATION bin)