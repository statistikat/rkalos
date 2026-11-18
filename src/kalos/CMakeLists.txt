cmake_minimum_required(VERSION 3.20)
project(kalos VERSION 0.1.0 LANGUAGES CXX)
option(KALOS_OPENMP "Enable owned OpenMP regions when available" ON)
option(KALOS_BUILD_TESTS "Build standalone tests" ON)
option(KALOS_SANITIZERS "Enable address and undefined behavior sanitizers" OFF)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(KalosOpenBLAS REQUIRED)
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES Kalos::OpenBLAS)
check_cxx_source_compiles("#include <cblas.h>
#include <lapacke.h>
#include <type_traits>
static_assert(sizeof(int)==4 && std::is_same<lapack_int,int>::value, \"LP64 required\");
int main() { return openblas_get_config() ? 0 : 1; }" KALOS_LP64_PROVIDER)
unset(CMAKE_REQUIRED_LIBRARIES)
if(NOT KALOS_LP64_PROVIDER)
  message(FATAL_ERROR "Kalos requires OpenBLAS with CBLAS, LAPACKE and a 32-bit LP64 integer ABI")
endif()
# A committed binding snapshot supplies its own exact revision header.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/kalos/revision.hpp")
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/include/kalos/revision.hpp" KALOS_REVISION_HEADER)
else()
  if(NOT KALOS_REVISION)
    execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE KALOS_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  endif()
  if(NOT KALOS_REVISION)
    set(KALOS_REVISION "unrecorded-source-archive")
  endif()
  set(KALOS_REVISION_HEADER "#define KALOS_CORE_REVISION \"${KALOS_REVISION}\"\n")
endif()
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated/kalos")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/generated/kalos/revision.hpp" "${KALOS_REVISION_HEADER}")
add_library(kalos STATIC
  src/arma_identification.cpp
  src/arma_information.cpp
  src/arma_objective.cpp
  src/arma_optimizer.cpp
  src/arma_parameters.cpp
  src/arma_residuals.cpp
  src/differencing_decision.cpp
  src/diffuse_state_space.cpp
  src/estimation_control.cpp
  src/fdlibm.cpp
  src/finite_state_factor.cpp
  src/forecast_arma_information.cpp
  src/forecast_information.cpp
  src/hannan_rissanen.cpp
  src/information.cpp
  src/initial_arma_identification.cpp
  src/interpolation_arma_information.cpp
  src/interpolation_information.cpp
  src/joint_regression_information.cpp
  src/linear_algebra.cpp
  src/numerics.cpp
  src/optimizer.cpp
  src/outlier_scan_engine.cpp
  src/outlier_scoring.cpp
  src/prediction_arma_information.cpp
  src/preparation.cpp
  src/prediction_moments.cpp
  src/regression_likelihood.cpp
  src/residual_information.cpp
  src/sarima_stabilization.cpp
  src/sparse_history.cpp
  src/ssq_levmar.cpp
  src/ssq_deviances.cpp
  src/state_space.cpp
  src/transformation_decision.cpp
  src/transforms.cpp
  src/vector_jacobian.cpp
  src/x11_b_pass.cpp
  src/x11_c_pass.cpp
  src/x11_d_pass.cpp
  src/x11_filters.cpp
  src/x11_initial.cpp
  src/x11_mean_endpoints.cpp
  src/x11_msr.cpp
  src/x11_msr_selection.cpp
  src/x11_musgrave.cpp
  src/x11_normalization.cpp
  src/x11_selection_pass.cpp
  src/x11_seasonal_endpoints.cpp
  src/x11_seasonal_plan.cpp
  src/x11_seasonal_processing.cpp
  src/x11_replacement.cpp
  src/x11_sigma.cpp
  src/x11_trend_selection.cpp
  src/x13_allocation.cpp
  src/x13_fit.cpp
  backends/openblas.cpp)
add_library(kalos::kalos ALIAS kalos)
set_target_properties(kalos PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(kalos PUBLIC cxx_std_17)
target_include_directories(kalos PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(kalos PUBLIC Kalos::OpenBLAS)
set(KALOS_WITH_OPENMP OFF)
if(KALOS_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_LIBRARIES Kalos::OpenBLAS OpenMP::OpenMP_CXX)
    set(_kalos_probe_signature "${KalosOpenBLAS_LIBRARY};${OpenMP_CXX_FLAGS};${OpenMP_CXX_LIBRARIES}")
    if(NOT "${_kalos_probe_signature}" STREQUAL "${KALOS_OPENMP_PROBE_SIGNATURE}")
      unset(KALOS_OPENMP_RUNTIME_OK CACHE)
      set(KALOS_OPENMP_PROBE_SIGNATURE "${_kalos_probe_signature}" CACHE INTERNAL "OpenMP provider probe inputs" FORCE)
    endif()
    if(NOT CMAKE_CROSSCOMPILING)
      check_cxx_source_runs("#include <cblas.h>
#include <omp.h>
int main() {
  if (!openblas_get_config()) return 1;
  int count = 0;
  #pragma omp parallel num_threads(2) reduction(+:count)
  count += 1;
  return count < 1;
}" KALOS_OPENMP_RUNTIME_OK)
    endif()
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(KALOS_OPENMP_RUNTIME_OK)
      target_link_libraries(kalos PUBLIC OpenMP::OpenMP_CXX)
      set(KALOS_WITH_OPENMP ON)
    else()
      message(WARNING "OpenMP/provider runtime probe failed or was not run; using serial Kalos. Select one compatible runtime (for example OpenMP_omp_LIBRARY), or provide KALOS_OPENMP_RUNTIME_OK for a qualified cross toolchain.")
    endif()
  else()
    message(STATUS "OpenMP unavailable: using serial execution")
  endif()
endif()
if(MSVC)
  target_compile_options(kalos PRIVATE /W4)
else()
  target_compile_options(kalos PRIVATE -Wall -Wextra -Wpedantic)
endif()
if(KALOS_SANITIZERS)
  if(MSVC)
    message(FATAL_ERROR "KALOS_SANITIZERS requires a Clang or GCC toolchain")
  endif()
  target_compile_options(kalos PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(kalos PUBLIC -fsanitize=address,undefined)
endif()
if(KALOS_BUILD_TESTS)
  enable_testing()
  foreach(test_name IN ITEMS
    test_arma_identification
    test_arma_information
    test_arma_objective
    test_arma_optimizer
    test_arma_parameters
    test_arma_residuals
    test_differencing_decision
    test_diffuse_state_space
    test_estimation_control
    test_forecast_arma_information
    test_forecast_information
    test_hannan_rissanen
    test_information
    test_initial_arma_identification
    test_interpolation_arma_information
    test_interpolation_information
    test_joint_regression_information
    test_legacy_arma_residuals
    test_linear_algebra
    test_numerics
    test_optimizer
    test_outlier_scoring
    test_parallel
    test_preparation
    test_prediction_moments
    test_regression_information
    test_regression_likelihood
    test_residual_information
    test_state_space
    test_transformation_decision
    test_transforms
    test_vector_jacobian
    test_x11_b_pass
    test_x11_c_pass
    test_x11_d_pass
    test_x11_filters
    test_x11_initial
    test_x11_log_additive
    test_x11_mean_endpoints
    test_x11_msr
    test_x11_msr_selection
    test_x11_multiplicative
    test_x11_musgrave
    test_x11_normalization
    test_x11_pseudo_additive
    test_x11_selection_pass
    test_x11_seasonal_endpoints
    test_x11_seasonal_plan
    test_x11_seasonal_processing
    test_x11_replacement
    test_x11_sigma
    test_x11_trend_selection
    test_x13_allocation
    test_x13_fit)
    add_executable(${test_name} tests/${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE kalos::kalos)
    # These inherited harnesses use assert; keep assertions in Release builds.
    if(MSVC)
      target_compile_options(${test_name} PRIVATE /UNDEBUG)
    else()
      target_compile_options(${test_name} PRIVATE -UNDEBUG)
    endif()
    add_test(NAME ${test_name} COMMAND ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT
      "KALOS_NUM_THREADS=2;OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
  endforeach()
endif()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(TARGETS kalos EXPORT kalosTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/kalos/revision.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/kalos)
install(EXPORT kalosTargets NAMESPACE kalos:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/kalos)
configure_package_config_file(cmake/kalosConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/kalosConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/kalos)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/kalosConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kalosConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/kalosConfigVersion.cmake cmake/FindKalosOpenBLAS.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/kalos)
install(FILES LICENSE.md DESTINATION ${CMAKE_INSTALL_DATADIR}/kalos)
