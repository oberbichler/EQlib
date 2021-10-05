
include(FindPackageHandleStandardArgs)

find_path(MKL_INCLUDE_DIR
    NAMES
        mkl.h
    PATHS
        ${MKLROOT}
        ${MKLROOT}/include
)

find_library(MKL_LIBRARY
    NAMES
        mkl_rt
    PATHS
        ${MKLROOT}
        ${MKLROOT}/lib
)

find_library(MKL_BINARY
    NAMES
        mkl_rt.dll
    PATHS
        ${MKLROOT}
        ${MKLROOT}/bin
)

find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_LIBRARY MKL_INCLUDE_DIR
)

mark_as_advanced(MKL_LIBRARY MKL_INCLUDE_DIR)

if(MKL_FOUND AND NOT TARGET MKL::MKL)
    add_library(MKL::MKL SHARED IMPORTED)
    set_property(TARGET MKL::MKL PROPERTY IMPORTED_IMPLIB ${MKL_LIBRARY})
    target_include_directories(MKL::MKL INTERFACE ${MKL_INCLUDE_DIR})
endif()