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
set(_kalos_search_options "")
if(KalosOpenBLAS_ROOT)
  set(KalosOpenBLAS_PROVIDER_SOURCE "explicit-root")
  set(_kalos_include_hints "${KalosOpenBLAS_ROOT}")
  set(_kalos_library_hints "${KalosOpenBLAS_ROOT}")
  set(_kalos_search_options NO_DEFAULT_PATH)
else()
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(PC_OPENBLAS QUIET openblas)
  endif()
  set(_kalos_include_hints ${PC_OPENBLAS_INCLUDE_DIRS})
  set(_kalos_library_hints ${PC_OPENBLAS_LIBRARY_DIRS})
endif()

# Drop cache entries produced by the original find_path/find_library calls.
# Command-line overrides have their own cache help string.
get_property(_kalos_help CACHE KalosOpenBLAS_INCLUDE_DIR PROPERTY HELPSTRING)
if(_kalos_help STREQUAL "Path to a file.")
  unset(KalosOpenBLAS_INCLUDE_DIR CACHE)
endif()
get_property(_kalos_help CACHE KalosOpenBLAS_LIBRARY PROPERTY HELPSTRING)
if(_kalos_help STREQUAL "Path to a library.")
  unset(KalosOpenBLAS_LIBRARY CACHE)
endif()

# Keep explicit overrides separate from discovery. Re-search on each configure,
# including when a provider is replaced in place.
if(NOT KalosOpenBLAS_INCLUDE_DIR)
  unset(_kalos_headers CACHE)
  find_path(_kalos_headers lapacke.h HINTS ${_kalos_include_hints}
    ${_kalos_search_options} PATH_SUFFIXES include include/openblas openblas)
  set(KalosOpenBLAS_INCLUDE_DIR "${_kalos_headers}")
endif()
if(NOT KalosOpenBLAS_LIBRARY)
  unset(_kalos_library CACHE)
  find_library(_kalos_library NAMES openblas HINTS ${_kalos_library_hints}
    ${_kalos_search_options} PATH_SUFFIXES lib lib64)
  set(KalosOpenBLAS_LIBRARY "${_kalos_library}")
endif()
set(_kalos_paths_valid TRUE)
foreach(_path IN ITEMS "${KalosOpenBLAS_INCLUDE_DIR}/cblas.h"
    "${KalosOpenBLAS_INCLUDE_DIR}/lapacke.h" "${KalosOpenBLAS_LIBRARY}")
  if(NOT EXISTS "${_path}")
    set(_kalos_paths_valid FALSE)
  elseif(KalosOpenBLAS_ROOT)
    file(REAL_PATH "${_path}" _resolved)
    file(REAL_PATH "${KalosOpenBLAS_ROOT}" _root)
    cmake_path(IS_PREFIX _root "${_resolved}" NORMALIZE _inside)
    if(NOT _inside)
      set(_kalos_paths_valid FALSE)
    endif()
  endif()
endforeach()

unset(_kalos_lp64 CACHE)
function(_kalos_check_provider)
  # Toolchain files may select STATIC_LIBRARY for compile-only probes. This
  # check must also resolve the numerical symbols.
  set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)
  include(CMakePushCheckState)
  include(CheckCXXSourceCompiles)
  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_INCLUDES "${KalosOpenBLAS_INCLUDE_DIR}")
  set(CMAKE_REQUIRED_LIBRARIES "${KalosOpenBLAS_LIBRARY}")
  check_cxx_source_compiles([=[
#include <cblas.h>
#include <lapacke.h>
#include <type_traits>
static_assert(sizeof(int) == 4 && std::is_same<lapack_int,int>::value
              && sizeof(blasint) == 4, "LP64 required");
// Keep an actual link reference without executing a numerical kernel.
template<class T> void require_export(T address) {
  T volatile reference = address;
  (void)reference;
}
int main() {
  require_export(&cblas_dgemm);
  require_export(&LAPACKE_dgelsy_work);
  require_export(&LAPACKE_dgeev_work);
  require_export(&LAPACKE_dgeqrf_work);
  require_export(&LAPACKE_dgesvd_work);
  require_export(&LAPACKE_dpocon_work);
  require_export(&LAPACKE_dpotrf_work);
  require_export(&LAPACKE_dpotrs_work);
  require_export(&LAPACKE_dsyev_work);
  require_export(&LAPACKE_dlarfg_work);
  return openblas_get_config() == nullptr;
}
]=] _kalos_lp64)
  cmake_pop_check_state()
endfunction()
if(_kalos_paths_valid)
  _kalos_check_provider()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KalosOpenBLAS REQUIRED_VARS
  KalosOpenBLAS_INCLUDE_DIR KalosOpenBLAS_LIBRARY _kalos_paths_valid _kalos_lp64)
if(KalosOpenBLAS_FOUND AND NOT TARGET Kalos::OpenBLAS)
  add_library(Kalos::OpenBLAS UNKNOWN IMPORTED)
  set_target_properties(Kalos::OpenBLAS PROPERTIES
    IMPORTED_LOCATION "${KalosOpenBLAS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${KalosOpenBLAS_INCLUDE_DIR}")
endif()
