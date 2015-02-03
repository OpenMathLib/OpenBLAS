##
## Author: Hank Anderson <hank@statease.com>
## Copyright: (c) Stat-Ease, Inc.
## Created: 12/29/14
## Last Modified: 12/29/14
## Description: Ported from the OpenBLAS/f_check perl script.
##              This is triggered by prebuild.cmake and runs before any of the code is built.
##              Appends Fortran information to config.h and Makefile.conf.

# CMake vars set by this file:
# F_COMPILER
# FC
# BU
# NOFORTRAN
# NEED2UNDERSCORES
# FEXTRALIB

# Defines set by this file:
# BUNDERSCORE
# NEEDBUNDERSCORE
# NEED2UNDERSCORES

if (NOT ONLY_CBLAS)
  # N.B. f_check is not cross-platform, so instead try to use CMake variables
  # run f_check (appends to TARGET files)
#  message(STATUS "Running f_check...")
#  execute_process(COMMAND perl f_check ${TARGET_MAKE} ${TARGET_CONF} ${CMAKE_Fortran_COMPILER}
#    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

  # TODO: is BU makefile macro needed?
  # TODO: detect whether underscore needed, set #defines appropriately - use try_compile
  # TODO: set FEXTRALIB flags a la f_check?

  set(BU "_")
  file(APPEND ${TARGET_CONF}
    "#define BUNDERSCORE _\n"
    "#define NEEDBUNDERSCORE 1\n"
    "#define NEED2UNDERSCORES 0\n")

else ()

  #When we only build CBLAS, we set NOFORTRAN=2
  set(NOFORTRAN 2)
  set(NO_FBLAS 1)
  #set(F_COMPILER GFORTRAN) # CMake handles the fortran compiler
  set(BU "_")
  file(APPEND ${TARGET_CONF}
    "#define BUNDERSCORE _\n"
    "#define NEEDBUNDERSCORE 1\n")
endif()
