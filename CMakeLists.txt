##
## Author: Hank Anderson <hank@statease.com>
##

cmake_minimum_required(VERSION 2.8.5)
project(OpenBLAS)
set(OpenBLAS_MAJOR_VERSION 0)
set(OpenBLAS_MINOR_VERSION 2)
set(OpenBLAS_PATCH_VERSION 20)
set(OpenBLAS_VERSION "${OpenBLAS_MAJOR_VERSION}.${OpenBLAS_MINOR_VERSION}.${OpenBLAS_PATCH_VERSION}")

enable_language(ASM)
enable_language(C)

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

if(MSVC)
set(OpenBLAS_LIBNAME libopenblas)
else()
set(OpenBLAS_LIBNAME openblas)
endif()

#######
if(MSVC)
option(BUILD_WITHOUT_LAPACK "Without LAPACK and LAPACKE (Only BLAS or CBLAS)" ON)
endif()
option(BUILD_WITHOUT_CBLAS "Without CBLAS" OFF)
option(BUILD_DEBUG "Build Debug Version" OFF)
#######
if(BUILD_WITHOUT_LAPACK)
set(NO_LAPACK 1)
set(NO_LAPACKE 1)
endif()

if(CMAKE_CONFIGURATION_TYPES) # multiconfig generator?
        set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)
        set(CMAKE_BUILD_TYPE
	        Debug Debug
		Release Release
	)
else()
	if( NOT CMAKE_BUILD_TYPE )
		if(BUILD_DEBUG)
			set(CMAKE_BUILD_TYPE Debug)
		else()
			set(CMAKE_BUILD_TYPE Release)
		endif()
	endif()
endif()

if(BUILD_WITHOUT_CBLAS)
set(NO_CBLAS 1)
endif()

#######


message(WARNING "CMake support is experimental. This will not produce the same Makefiles that OpenBLAS ships with. Only x86 support is currently available.")

include("${PROJECT_SOURCE_DIR}/cmake/utils.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/system.cmake")

set(BLASDIRS interface driver/level2 driver/level3 driver/others)

if (NOT DYNAMIC_ARCH)
  list(APPEND BLASDIRS kernel)
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

#Set default output directory
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
if(MSVC)
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/lib/Debug)
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/lib/Release)
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
if (NOT NOFORTRAN AND NOT NO_LAPACK)
  include("${PROJECT_SOURCE_DIR}/cmake/lapack.cmake")
if (NOT NO_LAPACKE)
  include("${PROJECT_SOURCE_DIR}/cmake/lapacke.cmake")
endif ()
endif ()

# Only generate .def for dll on MSVC and always produce pdb files for debug and release
if(MSVC)
set(OpenBLAS_DEF_FILE "${PROJECT_BINARY_DIR}/openblas.def")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /Zi")
set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DEBUG /OPT:REF /OPT:ICF")
endif()

# add objects to the openblas lib
add_library(${OpenBLAS_LIBNAME} SHARED ${LA_SOURCES} ${LAPACKE_SOURCES} ${TARGET_OBJS} ${OpenBLAS_DEF_FILE})

include("${PROJECT_SOURCE_DIR}/cmake/export.cmake")

# Set output for libopenblas
set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES LIBRARY_OUTPUT_NAME_DEBUG "${OpenBLAS_LIBNAME}_d")

foreach (OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
  string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  
  set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib/${OUTPUTCONFIG} )
  set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib/${OUTPUTCONFIG} )
  set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib/${OUTPUTCONFIG} )
endforeach()

enable_testing()
add_subdirectory(utest)

if (NOT MSVC)
	#only build shared library for MSVC

	add_library(${OpenBLAS_LIBNAME}_static STATIC ${LA_SOURCES} ${LAPACKE_SOURCES} ${TARGET_OBJS})
	set_target_properties(${OpenBLAS_LIBNAME}_static PROPERTIES OUTPUT_NAME ${OpenBLAS_LIBNAME})
	set_target_properties(${OpenBLAS_LIBNAME}_static PROPERTIES CLEAN_DIRECT_OUTPUT 1)

	if(SMP)
	target_link_libraries(${OpenBLAS_LIBNAME} pthread)
	target_link_libraries(${OpenBLAS_LIBNAME}_static pthread)
endif()

#build test and ctest
add_subdirectory(test)
if(NOT NO_CBLAS)
add_subdirectory(ctest)
endif()
endif()

set_target_properties(${OpenBLAS_LIBNAME} PROPERTIES 
  VERSION ${OpenBLAS_MAJOR_VERSION}.${OpenBLAS_MINOR_VERSION}
  SOVERSION ${OpenBLAS_MAJOR_VERSION}
)


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

# Install project

# Install libraries
install(TARGETS ${OpenBLAS_LIBNAME}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} )

# Install include files
	set (GENCONFIG_BIN ${CMAKE_BINARY_DIR}/gen_config_h${CMAKE_EXECUTABLE_SUFFIX})
	ADD_CUSTOM_COMMAND(
	OUTPUT ${CMAKE_BINARY_DIR}/openblas_config.h
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/config.h
	COMMAND ${GENCONFIG_BIN} ${CMAKE_CURRENT_SOURCE_DIR}/config.h ${CMAKE_CURRENT_SOURCE_DIR}/openblas_config_template.h > ${CMAKE_BINARY_DIR}/openblas_config.h
	)

	ADD_CUSTOM_TARGET(genconfig
	ALL
	DEPENDS openblas_config.h
	)
	add_dependencies(genconfig ${OpenBLAS_LIBNAME})

	install (FILES ${CMAKE_BINARY_DIR}/openblas_config.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

	message(STATUS "Generating f77blas.h in ${CMAKE_INSTALL_INCLUDEDIR}")

	ADD_CUSTOM_TARGET(genf77blas
	ALL
	COMMAND ${AWK} 'BEGIN{print \"\#ifndef OPENBLAS_F77BLAS_H\" \; print \"\#define OPENBLAS_F77BLAS_H\" \; print \"\#include \\"openblas_config.h\\" \"}; NF {print}; END{print \"\#endif\"}' ${CMAKE_CURRENT_SOURCE_DIR}/common_interface.h > ${CMAKE_BINARY_DIR}/f77blas.h
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/config.h
	)
	add_dependencies(genf77blas ${OpenBLAS_LIBNAME})

	install (FILES ${CMAKE_BINARY_DIR}/f77blas.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(NOT NO_CBLAS)
	message (STATUS "Generating cblas.h in ${CMAKE_INSTALL_INCLUDEDIR}")

	ADD_CUSTOM_TARGET(gencblas
	ALL
	COMMAND ${SED} 's/common/openblas_config/g' ${CMAKE_CURRENT_SOURCE_DIR}/cblas.h > "${CMAKE_BINARY_DIR}/cblas.tmp"
	COMMAND cp "${CMAKE_BINARY_DIR}/cblas.tmp" "${CMAKE_BINARY_DIR}/cblas.h"
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cblas.h
	)
	add_dependencies(gencblas ${OpenBLAS_LIBNAME})

	install (FILES ${CMAKE_BINARY_DIR}/cblas.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

if(NOT NO_LAPACKE)
	message (STATUS "Copying LAPACKE header files to ${CMAKE_INSTALL_INCLUDEDIR}")
	add_dependencies( ${OpenBLAS_LIBNAME} genlapacke)
	FILE(GLOB_RECURSE INCLUDE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lapack-netlib/LAPACKE/*.h")
	install (FILES ${INCLUDE_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
 
	ADD_CUSTOM_TARGET(genlapacke
	COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/lapack-netlib/LAPACKE/include/lapacke_mangling_with_flags.h.in "${CMAKE_BINARY_DIR}/lapacke_mangling.h"
	)
	install (FILES ${CMAKE_BINARY_DIR}/lapacke_mangling.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()
 
if(NOT MSVC)
	install (TARGETS ${OpenBLAS_LIBNAME}_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

include(FindPkgConfig QUIET)
if(PKG_CONFIG_FOUND)
	configure_file(${PROJECT_SOURCE_DIR}/cmake/openblas.pc.in ${PROJECT_BINARY_DIR}/openblas.pc @ONLY)
	install (FILES ${PROJECT_BINARY_DIR}/openblas.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/)
endif()
