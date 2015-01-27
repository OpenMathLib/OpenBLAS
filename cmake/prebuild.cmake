##
## Author: Hank Anderson <hank@statease.com>
## Copyright: (c) Stat-Ease, Inc.
## Created: 12/29/14
## Last Modified: 12/29/14
## Description: Ported from OpenBLAS/Makefile.prebuild
##              This is triggered by system.cmake and runs before any of the code is built.
##              Creates config.h and Makefile.conf by first running the c_check perl script (which creates those files).
##              Next it runs f_check and appends some fortran information to the files.
##              Finally it runs getarch and getarch_2nd for even more environment information.

# List of vars set by this file and included files:
# OSNAME
# ARCH
# C_COMPILER
# BINARY32
# BINARY64
# CEXTRALIB
# F_COMPILER
# FC
# BU
# CORE <- REQUIRED
# LIBCORE
# NUM_CORES <- REQUIRED
# HAVE_MMX
# HAVE_SSE
# HAVE_SSE2
# HAVE_SSE3
# MAKE
# SGEMM_UNROLL_M
# SGEMM_UNROLL_N
# DGEMM_UNROLL_M
# DGEMM_UNROLL_M
# QGEMM_UNROLL_N
# QGEMM_UNROLL_N
# CGEMM_UNROLL_M
# CGEMM_UNROLL_M
# ZGEMM_UNROLL_N
# ZGEMM_UNROLL_N
# XGEMM_UNROLL_M
# XGEMM_UNROLL_N
# CGEMM3M_UNROLL_M
# CGEMM3M_UNROLL_N
# ZGEMM3M_UNROLL_M
# ZGEMM3M_UNROLL_M
# XGEMM3M_UNROLL_N
# XGEMM3M_UNROLL_N

# CPUIDEMU = ../../cpuid/table.o

if (DEFINED CPUIDEMU)
  set(EXFLAGS "-DCPUIDEMU -DVENDOR=99")
endif ()

if (DEFINED TARGET_CORE)
  # set the C flags for just this file
  set_source_files_properties(getarch_2nd.c PROPERTIES COMPILE_FLAGS "-DBUILD_KERNEL")
  set(TARGET_MAKE "Makefile_kernel.conf")
  set(TARGET_CONF "config_kernel.h")
else()
  set(TARGET_MAKE "Makefile.conf")
  set(TARGET_CONF "config.h")
endif ()

include("${CMAKE_SOURCE_DIR}/cmake/c_check.cmake")
include("${CMAKE_SOURCE_DIR}/cmake/f_check.cmake")

# compile getarch
enable_language(ASM)
set(GETARCH_DIR "${PROJECT_BINARY_DIR}/getarch_build")
file(MAKE_DIRECTORY ${GETARCH_DIR})
try_compile(GETARCH_RESULT ${GETARCH_DIR}
  SOURCES ${CMAKE_SOURCE_DIR}/getarch.c ${CMAKE_SOURCE_DIR}/cpuid.S ${CPUIDEMO}
  COMPILE_DEFINITIONS ${EXFLAGS} -I${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GETARCH_LOG
  )

message(STATUS "GETARCH RESULT: ${GETARCH_RESULT}")
message(STATUS "GETARCH LOG: ${GETARCH_LOG}")

# TODO: need to append output of getarch binary to TARGET_CONF, not sure if I can get at it after using try_compile - may need to create CMakeLists.txt on the fly and build/execute

#add_executable(getarch getarch.c cpuid.S ${CPUIDEMU}
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
#
## run getarch, which appends even more to the TARGET files
#message(STATUS "Running getarch")
#execute_process(COMMAND getarch 0 >> ${TARGET_MAKE}
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
#execute_process(COMMAND getarch 1 >> ${TARGET_CONF}
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
#
## config.h is ready for getarch_2nd now, so compile that
#set(GETARCH2_SOURCES getarch_2nd.c config.h)
#add_executable(getarch_2nd getarch_2nd.c config.h)
#
## finally run getarch_2nd, appending yet more to the TARGET files
#message(STATUS "Running getarch_2nd")
#execute_process(COMMAND getarch_2nd 0 >> ${TARGET_MAKE}
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
#execute_process(COMMAND getarch_2nd 1 >> ${TARGET_CONF}
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

# TODO: need to read in the vars from Makefile.conf/Makefile_kernel.conf

# temporarily hardcoded to get system.cmake working
set(NUM_CORES 4)
set(CORE "GENERIC")

