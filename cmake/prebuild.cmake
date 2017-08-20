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

if(CMAKE_CROSSCOMPILING AND NOT DEFINED TARGET_CORE)
  # Detect target without running getarch
  if(AARCH64)
    set(TARGET_CORE "ARMV8")
  else()
    message(FATAL_ERROR "When cross compiling, a TARGET_CORE is required.")
  endif()
endif()

if (DEFINED TARGET_CORE)
  # set the C flags for just this file
  set(GETARCH2_FLAGS "-DBUILD_KERNEL")
  set(TARGET_MAKE "Makefile_kernel.conf")
  set(TARGET_CONF "config_kernel.h")
  set(TARGET_CONF_DIR ${PROJECT_BINARY_DIR}/kernel_config/${TARGET_CORE})
else()
  set(TARGET_MAKE "Makefile.conf")
  set(TARGET_CONF "config.h")
  set(TARGET_CONF_DIR ${PROJECT_BINARY_DIR})
endif ()

set(TARGET_CONF_TEMP "${PROJECT_BINARY_DIR}/${TARGET_CONF}.tmp")
include("${PROJECT_SOURCE_DIR}/cmake/c_check.cmake")

if (NOT NOFORTRAN)
  include("${PROJECT_SOURCE_DIR}/cmake/f_check.cmake")
endif ()

# Cannot run getarch on target if we are cross-compiling
if (CMAKE_CROSSCOMPILING)
  # Write to config as getarch would

  # TODO: Set up defines that getarch sets up based on every other target
  # Perhaps this should be inside a different file as it grows larger
  file(APPEND ${TARGET_CONF_TEMP}
    "#define ${TARGET_CORE}\n"
    "#define CHAR_CORENAME \"${TARGET_CORE}\"\n")
  if ("${TARGET_CORE}" STREQUAL "ARMV8")
    file(APPEND ${TARGET_CONF_TEMP}
      "#define L1_DATA_SIZE\t32768\n"
      "#define L1_DATA_LINESIZE\t64\n"
      "#define L2_SIZE\t262144\n"
      "#define L2_LINESIZE\t64\n"
      "#define DTB_DEFAULT_ENTRIES\t64\n"
      "#define DTB_SIZE\t4096\n"
      "#define L2_ASSOCIATIVE\t32\n")
    set(SGEMM_UNROLL_M 4)
    set(SGEMM_UNROLL_N 4)
  elseif ("${TARGET_CORE}" STREQUAL "CORTEXA57")
    file(APPEND ${TARGET_CONF_TEMP}
      "#define L1_CODE_SIZE\t49152\n"
      "#define L1_CODE_LINESIZE\t64\n"
      "#define L1_CODE_ASSOCIATIVE\t3\n"
      "#define L1_DATA_SIZE\t32768\n"
      "#define L1_DATA_LINESIZE\t64\n"
      "#define L1_DATA_ASSOCIATIVE\t2\n"
      "#define L2_SIZE\t2097152\n"
      "#define L2_LINESIZE\t64\n"
      "#define L2_ASSOCIATIVE\t16\n"
      "#define DTB_DEFAULT_ENTRIES\t64\n"
      "#define DTB_SIZE\t4096\n"
      "#define HAVE_VFPV4\n"
      "#define HAVE_VFPV3\n"
      "#define HAVE_VFP\n"
      "#define HAVE_NEON\n")
    set(SGEMM_DEFAULT_UNROLL_M 16)
    set(SGEMM_DEFAULT_UNROLL_N 4)
    set(DGEMM_DEFAULT_UNROLL_M 8)
    set(DGEMM_DEFAULT_UNROLL_N 4)
    set(CGEMM_DEFAULT_UNROLL_M 8)
    set(CGEMM_DEFAULT_UNROLL_N 4)
    set(ZGEMM_DEFAULT_UNROLL_M 8)
    set(ZGEMM_DEFAULT_UNROLL_N 4)
  endif()

  # Or should this actually be NUM_CORES?
  if (${NUM_THREADS} GREATER 0)
    file(APPEND ${TARGET_CONF_TEMP} "#define NUM_CORES\t${NUM_THREADS}\n")
  endif()

  # GetArch_2nd
  foreach(float_char S;D;Q;C;Z;X)
    if (NOT DEFINED ${float_char}GEMM_UNROLL_M)
      set(${float_char}GEMM_UNROLL_M 2)
    endif()
    if (NOT DEFINED ${float_char}GEMM_UNROLL_N)
      set(${float_char}GEMM_UNROLL_N 2)
    endif()
  endforeach()
  file(APPEND ${TARGET_CONF_TEMP}
    "#define GEMM_MULTITHREAD_THRESHOLD\t${GEMM_MULTITHREAD_THRESHOLD}\n")
  # Move to where gen_config_h would place it
  file(RENAME ${TARGET_CONF_TEMP} "${PROJECT_BINARY_DIR}/config.h")  

else()
# compile getarch
set(GETARCH_SRC
  ${PROJECT_SOURCE_DIR}/getarch.c
  ${CPUIDEMO}
)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  #Use generic for MSVC now
  message("MSVC")
  set(GETARCH_FLAGS ${GETARCH_FLAGS} -DFORCE_GENERIC)
else()
  list(APPEND GETARCH_SRC ${PROJECT_SOURCE_DIR}/cpuid.S)
endif ()

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "WindowsStore")
  # disable WindowsStore strict CRT checks
  set(GETARCH_FLAGS ${GETARCH_FLAGS} -D_CRT_SECURE_NO_WARNINGS)
endif ()

set(GETARCH_DIR "${PROJECT_BINARY_DIR}/getarch_build")
set(GETARCH_BIN "getarch${CMAKE_EXECUTABLE_SUFFIX}")
file(MAKE_DIRECTORY ${GETARCH_DIR})
configure_file(${TARGET_CONF_TEMP} ${GETARCH_DIR}/${TARGET_CONF} COPYONLY)
if (NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "WindowsStore")
  try_compile(GETARCH_RESULT ${GETARCH_DIR}
    SOURCES ${GETARCH_SRC}
    COMPILE_DEFINITIONS ${EXFLAGS} ${GETARCH_FLAGS} -I${GETARCH_DIR} -I${PROJECT_SOURCE_DIR} -I${PROJECT_BINARY_DIR}
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
file(APPEND ${TARGET_CONF_TEMP} ${GETARCH_CONF_OUT})
ParseGetArchVars(${GETARCH_MAKE_OUT})

set(GETARCH2_DIR "${PROJECT_BINARY_DIR}/getarch2_build")
set(GETARCH2_BIN "getarch_2nd${CMAKE_EXECUTABLE_SUFFIX}")
file(MAKE_DIRECTORY ${GETARCH2_DIR})
configure_file(${TARGET_CONF_TEMP} ${GETARCH2_DIR}/${TARGET_CONF} COPYONLY)
if (NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "WindowsStore")
  try_compile(GETARCH2_RESULT ${GETARCH2_DIR}
    SOURCES ${PROJECT_SOURCE_DIR}/getarch_2nd.c
    COMPILE_DEFINITIONS ${EXFLAGS} ${GETARCH_FLAGS} ${GETARCH2_FLAGS} -I${GETARCH2_DIR} -I${PROJECT_SOURCE_DIR} -I${PROJECT_BINARY_DIR}
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
file(APPEND ${TARGET_CONF_TEMP} ${GETARCH2_CONF_OUT})

if (${BUILD_KERNEL})
    configure_file(${TARGET_CONF_TEMP} ${PROJECT_BINARY_DIR}/kernel_config/${TARGET_CORE}/${TARGET_CONF} COPYONLY)
else ()
    configure_file(${TARGET_CONF_TEMP} ${PROJECT_BINARY_DIR}/${TARGET_CONF} COPYONLY)
endif ()

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

endif(CMAKE_CROSSCOMPILING)
