cmake_minimum_required(VERSION 3.16)

project(spip LANGUAGES C CXX)

option(SPIP_PREDICATES_USE_FLOAT "Use float version of Shewchuk predicates" OFF)

option(SPIP_BUILD_TESTS "Build spip tests" ON)

if(SPIP_BUILD_TESTS)
  enable_testing()
  add_executable(spip_test_pip
    tests/test_pip_basic.cpp
  )
  target_link_libraries(spip_test_pip PRIVATE spip)
  add_test(NAME spip_test_pip COMMAND spip_test_pip)

  add_executable(spip_test_psm_compile
  tests/test_psm_compile.cpp
  )
  target_link_libraries(spip_test_psm_compile PRIVATE spip)
  add_test(NAME spip_test_psm_compile COMMAND spip_test_psm_compile)

  add_executable(spip_test_quadruple3d
  tests/test_quadruple3d_basic.cpp
  )
  target_link_libraries(spip_test_quadruple3d PRIVATE spip)
  add_test(NAME spip_test_quadruple3d COMMAND spip_test_quadruple3d)
endif()

add_library(spip STATIC
  third_party/predicates.c
  third_party/geogram_psm/MultiPrecision_psm.cpp
  src/algorithms/point_in_polygon_sphere.cpp
  src/predicates/quadruple3d_exact.cpp
)

target_compile_features(spip PUBLIC cxx_std_17)

# Public headers live in include/, and orient3d.hpp includes "predicates.h"
# so we expose third_party/ as an include dir for consumers too.
target_include_directories(spip
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/geogram_psm>
    $<INSTALL_INTERFACE:include>
)

if(SPIP_PREDICATES_USE_FLOAT)
  target_compile_definitions(spip PUBLIC SPIP_PREDICATES_USE_FLOAT)
endif()

# Reasonable warnings
if(MSVC)
  target_compile_options(spip PRIVATE /W4)
else()
  target_compile_options(spip PRIVATE -Wall -Wextra -Wpedantic)
endif()
