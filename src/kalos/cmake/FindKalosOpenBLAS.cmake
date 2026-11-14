find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_OPENBLAS QUIET openblas)
endif()
find_path(KalosOpenBLAS_INCLUDE_DIR lapacke.h
  HINTS ${PC_OPENBLAS_INCLUDE_DIRS} ${OpenBLAS_ROOT}/include ENV OpenBLAS_ROOT
  PATH_SUFFIXES include include/openblas openblas)
find_library(KalosOpenBLAS_LIBRARY NAMES openblas
  HINTS ${PC_OPENBLAS_LIBRARY_DIRS} ${OpenBLAS_ROOT}/lib
  PATH_SUFFIXES lib lib64)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KalosOpenBLAS REQUIRED_VARS
  KalosOpenBLAS_INCLUDE_DIR KalosOpenBLAS_LIBRARY)
if(KalosOpenBLAS_FOUND AND NOT TARGET Kalos::OpenBLAS)
  add_library(Kalos::OpenBLAS UNKNOWN IMPORTED)
  set_target_properties(Kalos::OpenBLAS PROPERTIES
    IMPORTED_LOCATION "${KalosOpenBLAS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${KalosOpenBLAS_INCLUDE_DIR}")
endif()
