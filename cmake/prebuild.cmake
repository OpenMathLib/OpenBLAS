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
# TODO: need to use execute_process here, or compilation won't happen until later - maybe make temporary CMakeLists.txt file using file() ?
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

