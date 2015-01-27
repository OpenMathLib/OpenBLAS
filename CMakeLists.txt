##
## Author: Hank Anderson <hank@statease.com>
## Copyright: (c) Stat-Ease, Inc.
## Created: 12/23/14
## Last Modified: 12/23/14
##

cmake_minimum_required(VERSION 2.8.4)
project(OpenBLAS)

# is this necessary? lapack-netlib has its own fortran checks in its CMakeLists.txt
#enable_language(Fortran)

message(WARNING "CMake support is experimental. This will not produce the same Makefiles that OpenBLAS ships with.")

include("${CMAKE_SOURCE_DIR}/cmake/system.cmake")

set(BLASDIRS interface driver/level2 driver/level3 driver/others)

if (NOT ${DYNAMIC_ARCH})
  list(APPEND BLASDIRS kernel)
endif ()

if (DEFINED UTEST_CHECK)
  set(SANITY_CHECK 1)
endif ()

if (DEFINED SANITY_CHECK)
  list(APPEND BLASDIRS reference)
endif ()

set(SUBDIRS	${BLASDIRS})
if (NOT ${NO_LAPACK})
  list(APPEND SUBDIRS lapack)
endif ()

set(SUBDIRS_ALL ${SUBDIRS} test ctest utest exports benchmark ../laswp ../bench)

# all :: libs netlib tests shared

# libs:
if (NOT DEFINED CORE OR "${CORE}" STREQUAL "UNKNOWN")
  message(FATAL_ERROR "Detecting CPU failed. Please set TARGET explicitly, e.g. make TARGET=your_cpu_target. Please read README for details.")
endif ()

# Let CMake handle this
#if (${NOFORTRAN})
#  message(ERROR "OpenBLAS: Detecting fortran compiler failed. Please install fortran compiler, e.g. gfortran, ifort, openf90.")
#endif ()

if (${NO_STATIC} AND ${NO_SHARED})
  message(FATAL_ERROR "Neither static nor shared are enabled.")
endif ()

foreach (BLAS_DIR ${BLASDIRS})
  add_subdirectory(${BLAS_DIR})
endforeach ()

