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
  set(GETARCH2_FLAGS "-DBUILD_KERNEL")
  set(TARGET_MAKE "Makefile_kernel.conf")
  set(TARGET_CONF "config_kernel.h")
else()
  set(TARGET_MAKE "Makefile.conf")
  set(TARGET_CONF "config.h")
endif ()

include("${CMAKE_SOURCE_DIR}/cmake/c_check.cmake")
include("${CMAKE_SOURCE_DIR}/cmake/f_check.cmake")

# Reads string from getarch into CMake vars. Format of getarch vars is VARNAME=VALUE
function(ParseGetArchVars GETARCH_IN)
  string(REGEX MATCHALL "[0-9_a-zA-Z]+=[0-9_a-zA-Z]+" GETARCH_RESULT_LIST "${GETARCH_IN}")
  foreach (GETARCH_LINE ${GETARCH_RESULT_LIST})
    # split the line into var and value, then assign the value to a CMake var
    string(REGEX MATCHALL "[0-9_a-zA-Z]+" SPLIT_VAR "${GETARCH_LINE}")
    list(GET SPLIT_VAR 0 VAR_NAME)
    list(GET SPLIT_VAR 1 VAR_VALUE)
    set(${VAR_NAME} ${VAR_VALUE} PARENT_SCOPE)
  endforeach ()
endfunction ()

# compile getarch
enable_language(ASM)
set(GETARCH_DIR "${PROJECT_BINARY_DIR}/getarch_build")
set(GETARCH_BIN "getarch${CMAKE_EXECUTABLE_SUFFIX}")
file(MAKE_DIRECTORY ${GETARCH_DIR})
try_compile(GETARCH_RESULT ${GETARCH_DIR}
  SOURCES ${CMAKE_SOURCE_DIR}/getarch.c ${CMAKE_SOURCE_DIR}/cpuid.S ${CPUIDEMO}
  COMPILE_DEFINITIONS ${EXFLAGS} ${GETARCH_FLAGS} -I${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GETARCH_LOG
  COPY_FILE ${GETARCH_BIN}
)

message(STATUS "Running getarch")

# use the cmake binary w/ the -E param to run a shell command in a cross-platform way
execute_process(COMMAND ${GETARCH_BIN} 0 OUTPUT_VARIABLE GETARCH_MAKE_OUT)
execute_process(COMMAND ${GETARCH_BIN} 1 OUTPUT_VARIABLE GETARCH_CONF_OUT)

#message(STATUS "GETARCH results:\n${GETARCH_MAKE_OUT}")

# append config data from getarch to the TARGET file and read in CMake vars
file(APPEND ${TARGET_CONF} ${GETARCH_CONF_OUT})
ParseGetArchVars(${GETARCH_MAKE_OUT})

set(GETARCH2_DIR "${PROJECT_BINARY_DIR}/getarch2_build")
set(GETARCH2_BIN "getarch_2nd${CMAKE_EXECUTABLE_SUFFIX}")
file(MAKE_DIRECTORY ${GETARCH2_DIR})
try_compile(GETARCH2_RESULT ${GETARCH2_DIR}
  SOURCES ${CMAKE_SOURCE_DIR}/getarch_2nd.c
  COMPILE_DEFINITIONS ${EXFLAGS} ${GETARCH_FLAGS} ${GETARCH2_FLAGS} -I${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GETARCH2_LOG
  COPY_FILE ${GETARCH2_BIN}
)

# use the cmake binary w/ the -E param to run a shell command in a cross-platform way
execute_process(COMMAND ${GETARCH2_BIN} 0 OUTPUT_VARIABLE GETARCH2_MAKE_OUT)
execute_process(COMMAND ${GETARCH2_BIN} 1 OUTPUT_VARIABLE GETARCH2_CONF_OUT)

# append config data from getarch_2nd to the TARGET file and read in CMake vars
file(APPEND ${TARGET_CONF} ${GETARCH2_CONF_OUT})
ParseGetArchVars(${GETARCH2_MAKE_OUT})

