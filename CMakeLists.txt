cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# ---- Project ----

project(
  eqlib
  VERSION 2.0.0
  LANGUAGES CXX
)

# ---- Include guards

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# --- Get CPM

include(cmake/CPM.cmake)

# --- Add dependencies

CPMAddPackage("gh:TheLartians/PackageProject.cmake@1.6.0")

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

# --- Create library

add_library(${PROJECT_NAME} INTERFACE)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

target_include_directories(
  ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                            $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

# --- Create an installable target

packageProject(
  NAME ${PROJECT_NAME}
  VERSION ${PROJECT_VERSION}
  NAMESPACE ${PROJECT_NAME}
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
  INCLUDE_DESTINATION include/${PROJECT_NAME}-${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)