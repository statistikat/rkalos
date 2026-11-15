# Kalos OpenBLAS provider finder.
#
# Precedence rule: an explicit provider root (cache entry -DOpenBLAS_ROOT=<dir>
# or the OpenBLAS_ROOT environment variable) is EXCLUSIVE. When it is set, the
# header and the library are searched only under that root, so a configuration
# can never mix pinned provider headers with a system library found through
# pkg-config (for example a system OpenBLAS without LAPACKE exports). If the
# explicit root lacks either component, this module fails loudly instead of
# silently falling back to system hints. Without an explicit root, pkg-config
# and system hints are used as before.

set(KalosOpenBLAS_ROOT "${OpenBLAS_ROOT}")
if(NOT KalosOpenBLAS_ROOT AND DEFINED ENV{OpenBLAS_ROOT})
  set(KalosOpenBLAS_ROOT "$ENV{OpenBLAS_ROOT}")
endif()

set(KalosOpenBLAS_PROVIDER_SOURCE "system")
if(KalosOpenBLAS_ROOT)
  set(KalosOpenBLAS_PROVIDER_SOURCE "explicit-root")
  find_path(KalosOpenBLAS_INCLUDE_DIR lapacke.h
    HINTS "${KalosOpenBLAS_ROOT}" NO_DEFAULT_PATH
    PATH_SUFFIXES include include/openblas openblas)
  find_library(KalosOpenBLAS_LIBRARY NAMES openblas
    HINTS "${KalosOpenBLAS_ROOT}" NO_DEFAULT_PATH
    PATH_SUFFIXES lib lib64)
  if(NOT KalosOpenBLAS_INCLUDE_DIR OR NOT KalosOpenBLAS_LIBRARY)
    message(WARNING
      "OpenBLAS_ROOT='${KalosOpenBLAS_ROOT}' does not contain both lapacke.h (under include/) and libopenblas (under lib/); system provider hints are deliberately NOT used as a fallback. Unset OpenBLAS_ROOT to use pkg-config/system discovery.")
  endif()
else()
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(PC_OPENBLAS QUIET openblas)
  endif()
  find_path(KalosOpenBLAS_INCLUDE_DIR lapacke.h
    HINTS ${PC_OPENBLAS_INCLUDE_DIRS}
    PATH_SUFFIXES include include/openblas openblas)
  find_library(KalosOpenBLAS_LIBRARY NAMES openblas
    HINTS ${PC_OPENBLAS_LIBRARY_DIRS}
    PATH_SUFFIXES lib lib64)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KalosOpenBLAS REQUIRED_VARS
  KalosOpenBLAS_INCLUDE_DIR KalosOpenBLAS_LIBRARY)
if(KalosOpenBLAS_FOUND AND NOT TARGET Kalos::OpenBLAS)
  add_library(Kalos::OpenBLAS UNKNOWN IMPORTED)
  set_target_properties(Kalos::OpenBLAS PROPERTIES
    IMPORTED_LOCATION "${KalosOpenBLAS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${KalosOpenBLAS_INCLUDE_DIR}")
endif()
