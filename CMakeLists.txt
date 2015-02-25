##
## Author: Hank Anderson <hank@statease.com>
##

cmake_minimum_required(VERSION 2.8.4)
project(OpenBLAS)
set(OpenBLAS_MAJOR_VERSION 0)
set(OpenBLAS_MINOR_VERSION 2)
set(OpenBLAS_PATCH_VERSION 13)
set(OpenBLAS_VERSION "${OpenBLAS_MAJOR_VERSION}.${OpenBLAS_MINOR_VERSION}.${OpenBLAS_PATCH_VERSION}")

enable_language(ASM)

message(WARNING "CMake support is experimental. This will not produce the same Makefiles that OpenBLAS ships with. Only x86 support is currently available.")

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

# set which float types we want to build for
if (NOT DEFINED BUILD_SINGLE AND NOT DEFINED BUILD_DOUBLE AND NOT DEFINED BUILD_COMPLEX AND NOT DEFINED BUILD_COMPLEX16)
  # if none are defined, build for all
  set(BUILD_SINGLE true)
  set(BUILD_DOUBLE true)
  set(BUILD_COMPLEX true)
  set(BUILD_COMPLEX16 true)
endif ()

set(FLOAT_TYPES "")
if (BUILD_SINGLE)
  message(STATUS "Building Single Precision")
  list(APPEND FLOAT_TYPES "SINGLE") # defines nothing
endif ()

if (BUILD_DOUBLE)
  message(STATUS "Building Double Precision")
  list(APPEND FLOAT_TYPES "DOUBLE") # defines DOUBLE
endif ()

if (BUILD_COMPLEX)
  message(STATUS "Building Complex Precision")
  list(APPEND FLOAT_TYPES "COMPLEX") # defines COMPLEX
endif ()

if (BUILD_COMPLEX16)
  message(STATUS "Building Double Complex Precision")
  list(APPEND FLOAT_TYPES "ZCOMPLEX") # defines COMPLEX and DOUBLE
endif ()

set(SUBDIRS_ALL ${SUBDIRS} test ctest utest exports benchmark ../laswp ../bench)

# all :: libs netlib tests shared

# libs :
if (NOT DEFINED CORE OR "${CORE}" STREQUAL "UNKNOWN")
  message(FATAL_ERROR "Detecting CPU failed. Please set TARGET explicitly, e.g. make TARGET=your_cpu_target. Please read README for details.")
endif ()

if (${NO_STATIC} AND ${NO_SHARED})
  message(FATAL_ERROR "Neither static nor shared are enabled.")
endif ()

# get obj vars into format that add_library likes: $<TARGET_OBJS:objlib> (see http://www.cmake.org/cmake/help/v3.0/command/add_library.html)
set(TARGET_OBJS "")
foreach (SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
  string(REPLACE "/" "_" subdir_obj ${SUBDIR})
  list(APPEND TARGET_OBJS "$<TARGET_OBJECTS:${subdir_obj}>")
endforeach ()

# netlib:

# Can't just use lapack-netlib's CMake files, since they are set up to search for BLAS, build and install a binary. We just want to build a couple of lib files out of lapack and lapacke.
# Not using add_subdirectory here because lapack-netlib already has its own CMakeLists.txt. Instead include a cmake script with the sources we want.
include("${CMAKE_SOURCE_DIR}/cmake/lapack.cmake")
set_source_files_properties(${LA_SOURCES} PROPERTIES COMPILE_FLAGS "${LAPACK_FFLAGS}")

# add objects to the openblas lib
add_library(openblas ${LA_SOURCES} ${TARGET_OBJS})

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

