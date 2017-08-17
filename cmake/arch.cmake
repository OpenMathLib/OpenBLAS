##
## Author: Hank Anderson <hank@statease.com>
## Description: Ported from portion of OpenBLAS/Makefile.system
##              Sets various variables based on architecture.

if (${ARCH} STREQUAL "x86" OR ${ARCH} STREQUAL "x86_64")

  if (${ARCH} STREQUAL "x86")
    if (NOT BINARY)
      set(NO_BINARY_MODE 1)
    endif ()
  endif ()

  if (NOT NO_EXPRECISION)
    if (${F_COMPILER} MATCHES "GFORTRAN")
      # N.B. I'm not sure if CMake differentiates between GCC and LSB -hpa
      if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_C_COMPILER_ID} STREQUAL "LSB")
        set(EXPRECISION	1)
        set(CCOMMON_OPT "${CCOMMON_OPT} -DEXPRECISION -m128bit-long-double")
        set(FCOMMON_OPT	"${FCOMMON_OPT} -m128bit-long-double")
      endif ()
      if (${CMAKE_C_COMPILER_ID} STREQUAL "Clang")
        set(EXPRECISION	1)
        set(CCOMMON_OPT "${CCOMMON_OPT} -DEXPRECISION")
        set(FCOMMON_OPT	"${FCOMMON_OPT} -m128bit-long-double")
      endif ()
    endif ()
  endif ()
endif ()

if (${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
  set(CCOMMON_OPT "${CCOMMON_OPT} -wd981")
endif ()

if (USE_OPENMP)
  # USE_SIMPLE_THREADED_LEVEL3 = 1
  # NO_AFFINITY = 1
  find_package(OpenMP)
  if (OpenMP_FOUND)
    set(CCOMMON_OPT "${CCOMMON_OPT} ${OpenMP_C_FLAGS} -DUSE_OPENMP")
    set(FCOMMON_OPT "${FCOMMON_OPT} ${OpenMP_Fortran_FLAGS}")
  elseif (UNIX)
    set(USE_OPENMP 0)
  endif()
endif ()


if (DYNAMIC_ARCH)
  if (${ARCH} STREQUAL "x86")
    set(DYNAMIC_CORE KATMAI COPPERMINE NORTHWOOD PRESCOTT BANIAS CORE2 PENRYN DUNNINGTON NEHALEM ATHLON OPTERON OPTERON_SSE3 BARCELONA BOBCAT ATOM NANO)
  endif ()

  if (${ARCH} STREQUAL "x86_64")
    set(DYNAMIC_CORE PRESCOTT CORE2 PENRYN DUNNINGTON NEHALEM OPTERON OPTERON_SSE3 BARCELONA BOBCAT ATOM NANO)
    if (NOT NO_AVX)
      set(DYNAMIC_CORE ${DYNAMIC_CORE} SANDYBRIDGE BULLDOZER PILEDRIVER STEAMROLLER EXCAVATOR)
    endif ()
    if (NOT NO_AVX2)
      set(DYNAMIC_CORE ${DYNAMIC_CORE} HASWELL ZEN)
    endif ()
  endif ()

  if (NOT DYNAMIC_CORE)
    unset(DYNAMIC_ARCH)
  endif ()
endif ()

if (${ARCH} STREQUAL "ia64")
  set(NO_BINARY_MODE 1)
  set(BINARY_DEFINED 1)

  if (${F_COMPILER} MATCHES "GFORTRAN")
    if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
      # EXPRECISION	= 1
      # CCOMMON_OPT	+= -DEXPRECISION
    endif ()
  endif ()
endif ()

if (${ARCH} STREQUAL "mips64")
  set(NO_BINARY_MODE 1)
endif ()

if (${ARCH} STREQUAL "alpha")
  set(NO_BINARY_MODE 1)
  set(BINARY_DEFINED 1)
endif ()

if (${ARCH} STREQUAL "arm")
  set(NO_BINARY_MODE 1)
  set(BINARY_DEFINED 1)
endif ()

if (${ARCH} STREQUAL "arm64")
  set(NO_BINARY_MODE 1)
  set(BINARY_DEFINED 1)
endif ()

