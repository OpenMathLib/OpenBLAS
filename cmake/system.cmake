##
## Author: Hank Anderson <hank@statease.com>
## Copyright: (c) Stat-Ease, Inc.
## Created: 12/29/14
## Last Modified: 12/29/14
## Description: Ported from OpenBLAS/Makefile.system
##

set(NETLIB_LAPACK_DIR "${CMAKE_SOURCE_DIR}/lapack-netlib")

# TODO: Makefile.system detects Darwin (mac) and switches to clang here -hpa
# http://stackoverflow.com/questions/714100/os-detecting-makefile

# TODO: Makefile.system sets HOSTCC = $(CC) here if not already set -hpa

# TARGET_CORE will override TARGET which is used in DYNAMIC_ARCH=1.
if (DEFINED TARGET_CORE)
  set(TARGET ${TARGET_CORE})
endif ()

# Force fallbacks for 32bit
if (DEFINED BINARY AND DEFINED TARGET AND BINARY EQUAL 32)
  message(STATUS "Compiling a ${BINARY}-bit binary.")
  set(NO_AVX 1)
  if (${TARGET} STREQUAL "HASWELL" OR ${TARGET} STREQUAL "SANDYBRIDGE")
    set(TARGET "NEHALEM")
  endif ()
  if (${TARGET} STREQUAL "BULLDOZER" OR ${TARGET} STREQUAL "PILEDRIVER")
    set(TARGET "BARCELONA")
  endif ()
endif ()

if (DEFINED TARGET)
  message(STATUS "Targetting the ${TARGET} architecture.")
  set(GETARCH_FLAGS "-DFORCE_${TARGET}")
endif ()

if (${INTERFACE64})
  message(STATUS "Using 64-bit integers.")
  set(GETARCH_FLAGS	"${GETARCH_FLAGS} -DUSE64BITINT")
endif ()

if (NOT DEFINED GEMM_MULTITHREAD_THRESHOLD)
  set(GEMM_MULTITHREAD_THRESHOLD 4)
endif ()
message(STATUS "GEMM multithread threshold set to ${GEMM_MULTITHREAD_THRESHOLD}.")
set(GETARCH_FLAGS	"${GETARCH_FLAGS} -DGEMM_MULTITHREAD_THRESHOLD=${GEMM_MULTITHREAD_THRESHOLD}")

if (${NO_AVX})
  message(STATUS "Disabling Advanced Vector Extensions (AVX).")
  set(GETARCH_FLAGS "${GETARCH_FLAGS} -DNO_AVX")
endif ()

if (${NO_AVX2})
  message(STATUS "Disabling Advanced Vector Extensions 2 (AVX2).")
  set(GETARCH_FLAGS "${GETARCH_FLAGS} -DNO_AVX2")
endif ()

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  set(GETARCH_FLAGS "${GETARCH_FLAGS} -g")
endif ()

# TODO: let CMake handle this? -hpa
#if (${QUIET_MAKE})
#  set(MAKE "${MAKE} -s")
#endif()

if (NOT DEFINED NO_PARALLEL_MAKE)
  set(NO_PARALLEL_MAKE 0)
endif ()
set(GETARCH_FLAGS	"${GETARCH_FLAGS} -DNO_PARALLEL_MAKE=${NO_PARALLEL_MAKE}")

if (CMAKE_CXX_COMPILER STREQUAL loongcc)
  set(GETARCH_FLAGS	"${GETARCH_FLAGS} -static")
endif ()

#if don't use Fortran, it will only compile CBLAS.
if (${ONLY_CBLAS})
  set(NO_LAPACK 1)
else ()
  set(ONLY_CBLAS 0)
endif ()

include("${CMAKE_SOURCE_DIR}/cmake/prebuild.cmake")

if (NOT DEFINED NUM_THREADS)
  # TODO: NUM_CORES comes from `getarch.c` or `cpuid_x86.c`. This is built and executed above in `Makefile.prebuild`, and the results are in `Makefile.conf` and `Makefile_kernel.conf`. -hpa
  set(NUM_THREADS ${NUM_CORES})
endif ()

if (${NUM_THREADS} EQUAL 1)
  # TODO: was "override USE_THREAD = 0", do we need "override" here? -hpa
  set(USE_THREAD 0)
endif ()

if (DEFINED USE_THREAD)
  if (NOT ${USE_THREAD})
    unset(SMP)
  else ()
    set(SMP 1)
  endif ()
else ()
  # N.B. this is NUM_THREAD in Makefile.system which is probably a bug -hpa
  if (${NUM_THREADS} EQUAL 1)
    unset(SMP)
  else ()
    set(SMP 1)
  endif ()
endif ()

if (${SMP})
  message("SMP enabled.")
endif ()

if (NOT DEFINED NEED_PIC)
  set(NEED_PIC 1)
endif ()

# TODO: I think CMake should be handling all this stuff -hpa
unset(ARFLAGS)
set(CPP "${COMPILER} -E")
set(AR "${CROSS_SUFFIX}ar")
set(AS "$(CROSS_SUFFIX)as")
set(LD "$(CROSS_SUFFIX)ld")
set(RANLIB "$(CROSS_SUFFIX)ranlib")
set(NM "$(CROSS_SUFFIX)nm")
set(DLLWRAP "$(CROSS_SUFFIX)dllwrap")
set(OBJCOPY "$(CROSS_SUFFIX)objcopy")
set(OBJCONV "$(CROSS_SUFFIX)objconv")

