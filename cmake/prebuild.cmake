##
## Author: Hank Anderson <hank@statease.com>
## Description: Ported from OpenBLAS/Makefile.prebuild
##              This is triggered by system.cmake and runs before any of the code is built.
##              Creates config.h and Makefile.conf by first running the c_check perl script (which creates those files).
##              Next it runs f_check and appends some fortran information to the files.
##              Then it runs getarch and getarch_2nd for even more environment information.
##              Finally it builds gen_config_h for use at build time to generate config.h.

# CMake vars set by this file:
# CORE
# LIBCORE
# NUM_CORES
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

include("${PROJECT_SOURCE_DIR}/cmake/c_check.cmake")

if (NOT NOFORTRAN)
  include("${PROJECT_SOURCE_DIR}/cmake/f_check.cmake")
endif ()

# compile getarch
set(GETARCH_SRC
  ${PROJECT_SOURCE_DIR}/getarch.c
  ${CPUIDEMO}
)

if (NOT MSVC)
  list(APPEND GETARCH_SRC ${PROJECT_SOURCE_DIR}/cpuid.S)
endif ()

if (MSVC)
#Use generic for MSVC now
set(GETARCH_FLAGS ${GETARCH_FLAGS} -DFORCE_GENERIC)
endif()

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "WindowsStore")
  # disable WindowsStore strict CRT checks
  set(GETARCH_FLAGS ${GETARCH_FLAGS} -D_CRT_SECURE_NO_WARNINGS)
endif ()

set(GETARCH_DIR "${PROJECT_BINARY_DIR}/getarch_build")
set(GETARCH_BIN "getarch${CMAKE_EXECUTABLE_SUFFIX}")
file(MAKE_DIRECTORY ${GETARCH_DIR})
if (NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "WindowsStore")
  try_compile(GETARCH_RESULT ${GETARCH_DIR}
    SOURCES ${GETARCH_SRC}
    COMPILE_DEFINITIONS ${EXFLAGS} ${GETARCH_FLAGS} -I${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GETARCH_LOG
    COPY_FILE ${PROJECT_BINARY_DIR}/${GETARCH_BIN}
  )

  if (NOT ${GETARCH_RESULT})
    MESSAGE(FATAL_ERROR "Compiling getarch failed ${GETARCH_LOG}")
  endif ()
endif ()
message(STATUS "Running getarch")

# use the cmake binary w/ the -E param to run a shell command in a cross-platform way
execute_process(COMMAND ${PROJECT_BINARY_DIR}/${GETARCH_BIN} 0 OUTPUT_VARIABLE GETARCH_MAKE_OUT)
execute_process(COMMAND ${PROJECT_BINARY_DIR}/${GETARCH_BIN} 1 OUTPUT_VARIABLE GETARCH_CONF_OUT)

message(STATUS "GETARCH results:\n${GETARCH_MAKE_OUT}")

# append config data from getarch to the TARGET file and read in CMake vars
file(APPEND ${TARGET_CONF} ${GETARCH_CONF_OUT})
ParseGetArchVars(${GETARCH_MAKE_OUT})

set(GETARCH2_DIR "${PROJECT_BINARY_DIR}/getarch2_build")
set(GETARCH2_BIN "getarch_2nd${CMAKE_EXECUTABLE_SUFFIX}")
file(MAKE_DIRECTORY ${GETARCH2_DIR})
if (NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "WindowsStore")
  try_compile(GETARCH2_RESULT ${GETARCH2_DIR}
    SOURCES ${PROJECT_SOURCE_DIR}/getarch_2nd.c
    COMPILE_DEFINITIONS ${EXFLAGS} ${GETARCH_FLAGS} ${GETARCH2_FLAGS} -I${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GETARCH2_LOG
    COPY_FILE ${PROJECT_BINARY_DIR}/${GETARCH2_BIN}
  )

  if (NOT ${GETARCH2_RESULT})
    MESSAGE(FATAL_ERROR "Compiling getarch_2nd failed ${GETARCH2_LOG}")
  endif ()
endif ()

# use the cmake binary w/ the -E param to run a shell command in a cross-platform way
execute_process(COMMAND ${PROJECT_BINARY_DIR}/${GETARCH2_BIN} 0 OUTPUT_VARIABLE GETARCH2_MAKE_OUT)
execute_process(COMMAND ${PROJECT_BINARY_DIR}/${GETARCH2_BIN} 1 OUTPUT_VARIABLE GETARCH2_CONF_OUT)

# append config data from getarch_2nd to the TARGET file and read in CMake vars
file(APPEND ${TARGET_CONF} ${GETARCH2_CONF_OUT})
ParseGetArchVars(${GETARCH2_MAKE_OUT})

# compile get_config_h
set(GEN_CONFIG_H_DIR "${PROJECT_BINARY_DIR}/genconfig_h_build")
set(GEN_CONFIG_H_BIN "gen_config_h${CMAKE_EXECUTABLE_SUFFIX}")
set(GEN_CONFIG_H_FLAGS "-DVERSION=\"${OpenBLAS_VERSION}\"")
file(MAKE_DIRECTORY ${GEN_CONFIG_H_DIR})

if (NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "WindowsStore")
  try_compile(GEN_CONFIG_H_RESULT ${GEN_CONFIG_H_DIR}
    SOURCES ${PROJECT_SOURCE_DIR}/gen_config_h.c
    COMPILE_DEFINITIONS ${EXFLAGS} ${GETARCH_FLAGS} ${GEN_CONFIG_H_FLAGS} -I${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GEN_CONFIG_H_LOG
    COPY_FILE ${PROJECT_BINARY_DIR}/${GEN_CONFIG_H_BIN}
  )

  if (NOT ${GEN_CONFIG_H_RESULT})
    MESSAGE(FATAL_ERROR "Compiling gen_config_h failed ${GEN_CONFIG_H_LOG}")
  endif ()
endif ()