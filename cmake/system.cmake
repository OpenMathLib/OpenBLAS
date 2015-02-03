##
## Author: Hank Anderson <hank@statease.com>
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
  message(STATUS "SMP enabled.")
endif ()

if (NOT DEFINED NEED_PIC)
  set(NEED_PIC 1)
endif ()

# TODO: I think CMake should be handling all this stuff -hpa
unset(ARFLAGS)
set(CPP "${COMPILER} -E")
set(AR "${CROSS_SUFFIX}ar")
set(AS "${CROSS_SUFFIX}as")
set(LD "${CROSS_SUFFIX}ld")
set(RANLIB "${CROSS_SUFFIX}ranlib")
set(NM "${CROSS_SUFFIX}nm")
set(DLLWRAP "${CROSS_SUFFIX}dllwrap")
set(OBJCOPY "${CROSS_SUFFIX}objcopy")
set(OBJCONV "${CROSS_SUFFIX}objconv")

# OS dependent settings
include("${CMAKE_SOURCE_DIR}/cmake/os.cmake")

# Architecture dependent settings
include("${CMAKE_SOURCE_DIR}/cmake/arch.cmake")

# C Compiler dependent settings
include("${CMAKE_SOURCE_DIR}/cmake/cc.cmake")

# Fortran Compiler dependent settings
include("${CMAKE_SOURCE_DIR}/cmake/fc.cmake")

if (BINARY64)
  if (INTERFACE64)
    # CCOMMON_OPT += -DUSE64BITINT
  endif ()
endif ()

if (NEED_PIC)
  if (${CMAKE_C_COMPILER} STREQUAL "IBM")
    set(CCOMMON_OPT "${CCOMMON_OPT} -qpic=large")
  else ()
    set(CCOMMON_OPT "${CCOMMON_OPT} -fPIC")
  endif ()

  if (${CMAKE_Fortran_COMPILER} STREQUAL "SUN")
    set(FCOMMON_OPT "${FCOMMON_OPT} -pic")
  else ()
    set(FCOMMON_OPT "${FCOMMON_OPT} -fPIC")
  endif ()
endif ()

if (DYNAMIC_ARCH)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DDYNAMIC_ARCH")
endif ()

if (NO_LAPACK)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_LAPACK")
  #Disable LAPACK C interface
  set(NO_LAPACKE 1)
endif ()

if (NO_LAPACKE)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_LAPACKE")
endif ()

if (NO_AVX)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_AVX")
endif ()

if (${ARCH} STREQUAL "x86")
  set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_AVX")
endif ()

if (NO_AVX2)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_AVX2")
endif ()

if (SMP)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DSMP_SERVER")

  if (${ARCH} STREQUAL "mips64")
    if (NOT ${CORE} STREQUAL "LOONGSON3B")
      set(USE_SIMPLE_THREADED_LEVEL3 1)
    endif ()
  endif ()

  if (USE_OPENMP)
    # USE_SIMPLE_THREADED_LEVEL3 = 1
    # NO_AFFINITY = 1
    set(CCOMMON_OPT "${CCOMMON_OPT} -DUSE_OPENMP")
  endif ()

  if (BIGNUMA)
    set(CCOMMON_OPT "${CCOMMON_OPT} -DBIGNUMA")
  endif ()

endif ()

if (NO_WARMUP)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_WARMUP")
endif ()

if (CONSISTENT_FPCSR)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DCONSISTENT_FPCSR")
endif ()

# Only for development
# set(CCOMMON_OPT "${CCOMMON_OPT} -DPARAMTEST")
# set(CCOMMON_OPT "${CCOMMON_OPT} -DPREFETCHTEST")
# set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_SWITCHING")
# set(USE_PAPI 1)

if (USE_PAPI)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DUSE_PAPI")
  set(EXTRALIB "${EXTRALIB} -lpapi -lperfctr")
endif ()

if (DYNAMIC_THREADS)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DDYNAMIC_THREADS")
endif ()

set(CCOMMON_OPT "${CCOMMON_OPT} -DMAX_CPU_NUMBER=${NUM_THREADS}")

if (USE_SIMPLE_THREADED_LEVEL3)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DUSE_SIMPLE_THREADED_LEVEL3")
endif ()

if (DEFINED LIBNAMESUFFIX)
  set(LIBPREFIX "libopenblas_${LIBNAMESUFFIX}")
else ()
  set(LIBPREFIX "libopenblas")
endif ()

if (NOT DEFINED SYMBOLPREFIX)
  set(SYMBOLPREFIX "")
endif ()

if (NOT DEFINED SYMBOLSUFFIX)
  set(SYMBOLSUFFIX "")
endif ()

set(KERNELDIR	"${CMAKE_SOURCE_DIR}/kernel/${ARCH}")

# TODO: nead to convert these Makefiles
# include ${CMAKE_SOURCE_DIR}/cmake/${ARCH}.cmake

# TODO: Need to figure out how to get $(*F) in cmake
set(CCOMMON_OPT "${CCOMMON_OPT} -DASMNAME=${FU}$(*F) -DASMFNAME=${FU}$(*F)${BU} -DNAME=$(*F)${BU} -DCNAME=$(*F) -DCHAR_NAME=\"$(*F)${BU}\" -DCHAR_CNAME=\"$(*F)\"")

if (${CORE} STREQUAL "PPC440")
  set(CCOMMON_OPT "${CCOMMON_OPT} -DALLOC_QALLOC")
endif ()

if (${CORE} STREQUAL "PPC440FP2")
  set(STATIC_ALLOCATION 1)
endif ()

if (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  set(NO_AFFINITY 1)
endif ()

if (NOT ${ARCH} STREQUAL "x86_64" AND NOT ${ARCH} STREQUAL "x86" AND NOT ${CORE} STREQUAL "LOONGSON3B")
  set(NO_AFFINITY 1)
endif ()

if (NO_AFFINITY)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DNO_AFFINITY")
endif ()

if (FUNCTION_PROFILE)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DFUNCTION_PROFILE")
endif ()

if (HUGETLB_ALLOCATION)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DALLOC_HUGETLB")
endif ()

if (DEFINED HUGETLBFILE_ALLOCATION)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DALLOC_HUGETLBFILE -DHUGETLB_FILE_NAME=${HUGETLBFILE_ALLOCATION})")
endif ()

if (STATIC_ALLOCATION)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DALLOC_STATIC")
endif ()

if (DEVICEDRIVER_ALLOCATION)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DALLOC_DEVICEDRIVER -DDEVICEDRIVER_NAME=\"/dev/mapper\"")
endif ()

if (MIXED_MEMORY_ALLOCATION)
  set(CCOMMON_OPT "${CCOMMON_OPT} -DMIXED_MEMORY_ALLOCATION")
endif ()

if (${CMAKE_SYSTEM_NAME} STREQUAL "SunOS")
  set(TAR	gtar)
  set(PATCH	gpatch)
  set(GREP ggrep)
else ()
  set(TAR tar)
  set(PATCH patch)
  set(GREP grep)
endif ()

if (NOT DEFINED MD5SUM)
  set(MD5SUM md5sum)
endif ()

set(AWK awk)

set(REVISION "-r${OpenBLAS_VERSION}")
set(MAJOR_VERSION ${OpenBLAS_MAJOR_VERSION})

if (DEBUG)
  set(COMMON_OPT "${COMMON_OPT} -g")
endif ()

if (NOT DEFINED COMMON_OPT)
  set(COMMON_OPT "-O2")
endif ()


