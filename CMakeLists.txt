##
## Author: Hank Anderson <hank@statease.com>
##

cmake_minimum_required(VERSION 2.8.4)
project(OpenBLAS)
set(OpenBLAS_MAJOR_VERSION 0)
set(OpenBLAS_MINOR_VERSION 2)
set(OpenBLAS_PATCH_VERSION 13)
set(OpenBLAS_VERSION "${OpenBLAS_MAJOR_VERSION}.${OpenBLAS_MINOR_VERSION}.${OpenBLAS_PATCH_VERSION}")

enable_language(Fortran)
enable_language(ASM)

message(WARNING "CMake support is experimental. This will not produce the same Makefiles that OpenBLAS ships with.")

include("${CMAKE_SOURCE_DIR}/cmake/utils.cmake")
include("${CMAKE_SOURCE_DIR}/cmake/system.cmake")

set(BLASDIRS interface driver/level2 driver/level3 driver/others)

if (NOT DYNAMIC_ARCH)
  list(APPEND BLASDIRS kernel)
endif ()

if (DEFINED UTEST_CHECK)
  set(SANITY_CHECK 1)
endif ()

if (DEFINED SANITY_CHECK)
  list(APPEND BLASDIRS reference)
endif ()

set(SUBDIRS	${BLASDIRS})
if (NOT NO_LAPACK)
  list(APPEND SUBDIRS lapack)
endif ()

set(SUBDIRS_ALL ${SUBDIRS} test ctest utest exports benchmark ../laswp ../bench)

# all :: libs netlib tests shared

# libs :
if (NOT DEFINED CORE OR "${CORE}" STREQUAL "UNKNOWN")
  message(FATAL_ERROR "Detecting CPU failed. Please set TARGET explicitly, e.g. make TARGET=your_cpu_target. Please read README for details.")
endif ()

# Let CMake handle this
#if (${NOFORTRAN})
#  message(ERROR "OpenBLAS: Detecting fortran compiler failed. Please install fortran compiler, e.g. gfortran, ifort, openf90.")
#endif ()

if (${NO_STATIC} AND ${NO_SHARED})
  message(FATAL_ERROR "Neither static nor shared are enabled.")
endif ()

set(DBLAS_OBJS "")
foreach (SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach ()

# get obj vars into format that add_library likes: $<TARGET_OBJS:objlib> (see http://www.cmake.org/cmake/help/v3.0/command/add_library.html)
set(TARGET_OBJS "")
foreach (DBLAS_OBJ ${DBLAS_OBJS})
  get_target_property(PREV_DEFS ${DBLAS_OBJ} COMPILE_DEFINITIONS)
  set_target_properties(${DBLAS_OBJ} PROPERTIES COMPILE_DEFINITIONS "${PREV_DEFS};DOUBLE")
  list(APPEND TARGET_OBJS "$<TARGET_OBJECTS:${DBLAS_OBJ}>")
endforeach ()

# netlib:

# Can't just use lapack-netlib's CMake files, since they are set up to search for BLAS, build and install a binary. We just want to build a couple of lib files out of lapack and lapacke.

# N.B. if a source is added or removed you MUST re-run the cmake command manually; make will not do it.
file(GLOB DLA_SOURCES "${NETLIB_LAPACK_DIR}/SRC/d*.f")

add_library(DLA_OBJ OBJECT ${DLA_SOURCES})
set_target_properties(${DLA_OBJ} PROPERTIES COMPILE_FLAGS ${LAPACK_FFLAGS})
list(APPEND TARGET_OBJS "$<TARGET_OBJECTS:DLA_OBJ>")

# add objects to the openblas lib
add_library(openblas ${TARGET_OBJS})

# TODO: Why is the config saved here? Is this necessary with CMake?
#Save the config files for installation
#	@cp Makefile.conf Makefile.conf_last
#	@cp config.h config_last.h
#ifdef QUAD_PRECISION
#	@echo "#define QUAD_PRECISION">> config_last.h
#endif
#ifeq ($(EXPRECISION), 1)
#	@echo "#define EXPRECISION">> config_last.h
#endif
###
#ifeq ($(DYNAMIC_ARCH), 1)
#	@$(MAKE) -C kernel commonlibs || exit 1
#	@for d in $(DYNAMIC_CORE) ; \
#	do  $(MAKE) GOTOBLAS_MAKEFILE= -C kernel TARGET_CORE=$$d kernel || exit 1 ;\
#	done
#	@echo DYNAMIC_ARCH=1 >> Makefile.conf_last
#endif
#ifdef USE_THREAD
#	@echo USE_THREAD=$(USE_THREAD) >>  Makefile.conf_last
#endif
#	@touch lib.grd

