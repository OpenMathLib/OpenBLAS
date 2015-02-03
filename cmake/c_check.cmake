##
## Author: Hank Anderson <hank@statease.com>
## Description: Ported from the OpenBLAS/c_check perl script.
##              This is triggered by prebuild.cmake and runs before any of the code is built.
##              Creates config.h and Makefile.conf.

# N.B. c_check (and ctest.c) is not cross-platform, so instead try to use CMake variables.

# TODO: detect NEED_FU
set(NEED_FU 1)

# Convert CMake vars into the format that OpenBLAS expects
string(TOUPPER ${CMAKE_SYSTEM_NAME} HOST_OS)
if (${HOST_OS} STREQUAL "WINDOWS")
  set(HOST_OS WINNT)
endif ()

# added by hpa - check size of void ptr to detect 64-bit compile
if (NOT DEFINED BINARY)
  set(BINARY 32)
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(BINARY 64)
  endif ()
endif ()

# CMake docs define these:
# CMAKE_SYSTEM_PROCESSOR - The name of the CPU CMake is building for.
# CMAKE_HOST_SYSTEM_PROCESSOR - The name of the CPU CMake is running on.
set(HOST_ARCH ${CMAKE_SYSTEM_PROCESSOR})
if (${HOST_ARCH} STREQUAL "AMD64")
  set(HOST_ARCH "X86_64")
endif ()

# If you are using a 32-bit compiler on a 64-bit system CMAKE_SYSTEM_PROCESSOR will be wrong
if (${HOST_ARCH} STREQUAL "X86_64" AND BINARY EQUAL 32)
  set(HOST_ARCH X86)
endif ()

set(COMPILER_ID ${CMAKE_CXX_COMPILER_ID})
if (${COMPILER_ID} STREQUAL "GNU")
  set(COMPILER_ID "GCC")
endif ()

file(WRITE ${TARGET_CONF}
  "#define OS_${HOST_OS}\t1\n"
  "#define ARCH_${HOST_ARCH}\t1\n"
  "#define C_${COMPILER_ID}\t1\n"
  "#define __${BINARY}BIT__\t1\n"
  "#define FUNDERSCORE\t${NEED_FU}\n")

