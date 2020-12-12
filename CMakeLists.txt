##
## Author: Hank Anderson <hank@statease.com>
##

cmake_minimum_required(VERSION 2.8.5)
project(OpenBLAS C ASM)
set(OpenBLAS_MAJOR_VERSION 0)
set(OpenBLAS_MINOR_VERSION 3)
set(OpenBLAS_PATCH_VERSION 13)
set(OpenBLAS_VERSION "${OpenBLAS_MAJOR_VERSION}.${OpenBLAS_MINOR_VERSION}.${OpenBLAS_PATCH_VERSION}")

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

include(CMakePackageConfigHelpers)


#######
if(MSVC)
option(BUILD_WITHOUT_LAPACK "Do not build LAPACK and LAPACKE (Only BLAS or CBLAS)" ON)
endif()
option(BUILD_WITHOUT_CBLAS "Do not build the C interface (CBLAS) to the BLAS functions" OFF)
option(DYNAMIC_ARCH "Include support for multiple CPU targets, with automatic selection at runtime (x86/x86_64, aarch64 or ppc only)" OFF)
option(DYNAMIC_OLDER "Include specific support for older x86 cpu models (Penryn,Dunnington,Atom,Nano,Opteron) with DYNAMIC_ARCH" OFF)
option(BUILD_RELAPACK "Build with ReLAPACK (recursive implementation of several LAPACK functions on top of standard LAPACK)" OFF)
option(USE_LOCKING "Use locks even in single-threaded builds to make them callable from multiple threads" OFF)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
option(NO_AFFINITY "Disable support for CPU affinity masks to avoid binding processes from e.g. R or numpy/scipy to a single core" ON)
else()
set(NO_AFFINITY 1)
endif()
option(CPP_THREAD_SAFETY_TEST "Run a massively parallel DGEMM test to confirm thread safety of the library (requires OpenMP and about 1.3GB of RAM)" OFF)
option(CPP_THREAD_SAFETY_GEMV "Run a massively parallel DGEMV test to confirm thread safety of the library (requires OpenMP)" OFF)

# Add a prefix or suffix to all exported symbol names in the shared library.
# Avoids conflicts with other BLAS libraries, especially when using
# 64 bit integer interfaces in OpenBLAS.

set(SYMBOLPREFIX "" CACHE STRING  "Add a prefix to all exported symbol names in the shared library to avoid conflicts with other BLAS libraries" )
set(SYMBOLSUFFIX "" CACHE STRING  "Add a suffix to all exported symbol names in the shared library, e.g. _64 for INTERFACE64 builds" )
#######
if(BUILD_WITHOUT_LAPACK)
set(NO_LAPACK 1)
set(NO_LAPACKE 1)
endif()

if(BUILD_WITHOUT_CBLAS)
set(NO_CBLAS 1)
endif()

#######

if(MSVC AND MSVC_STATIC_CRT)
    set(CompilerFlags
            CMAKE_CXX_FLAGS
            CMAKE_CXX_FLAGS_DEBUG
            CMAKE_CXX_FLAGS_RELEASE
            CMAKE_C_FLAGS
            CMAKE_C_FLAGS_DEBUG
            CMAKE_C_FLAGS_RELEASE
            )
    foreach(CompilerFlag ${CompilerFlags})
      string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
    endforeach()
endif()

message(WARNING "CMake support is experimental. It does not yet support all build options and may not produce the same Makefiles that OpenBLAS ships with.")

include("${PROJECT_SOURCE_DIR}/cmake/utils.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/system.cmake")

set(OpenBLAS_LIBNAME openblas${SUFFIX64_UNDERSCORE})

set(BLASDIRS interface driver/level2 driver/level3 driver/others)

if (NOT DYNAMIC_ARCH)
  list(APPEND BLASDIRS kernel)
endif ()

if (DEFINED SANITY_CHECK)
  list(APPEND BLASDIRS reference)
endif ()

set(SUBDIRS	${BLASDIRS})
if (NOT NO_LAPACK)
  if(BUILD_RELAPACK)
    list(APPEND SUBDIRS relapack/src)
  endif()
  list(APPEND SUBDIRS lapack)
endif ()

if (NOT DEFINED BUILD_BFLOAT16)
 set (BUILD_BFLOAT16 false)
endif ()
# set which float types we want to build for
if (NOT DEFINED BUILD_SINGLE AND NOT DEFINED BUILD_DOUBLE AND NOT DEFINED BUILD_COMPLEX AND NOT DEFINED BUILD_COMPLEX16)
  # if none are defined, build for all
#  set(BUILD_BFLOAT16 true)
  set(BUILD_SINGLE true)
  set(BUILD_DOUBLE true)
  set(BUILD_COMPLEX true)
  set(BUILD_COMPLEX16 true)
endif ()

if (NOT DEFINED BUILD_MATGEN)
  set(BUILD_MATGEN true)
endif()

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

if (BUILD_BFLOAT16)
  message(STATUS "Building Half Precision")
  list(APPEND FLOAT_TYPES "BFLOAT16") # defines nothing
endif ()

if (NOT DEFINED CORE OR "${CORE}" STREQUAL "UNKNOWN")
  message(FATAL_ERROR "Detecting CPU failed. Please set TARGET explicitly, e.g. make TARGET=your_cpu_target. Please read README for details.")
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
  if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_LESS 3.4)
    set(OpenBLAS_DEF_FILE "${PROJECT_BINARY_DIR}/openblas.def")
  endif()
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /Zi")
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DEBUG /OPT:REF /OPT:ICF")
endif()

if (${DYNAMIC_ARCH})
  add_subdirectory(kernel)
  foreach(TARGET_CORE ${DYNAMIC_CORE})
    message("${TARGET_CORE}")
    list(APPEND TARGET_OBJS "$<TARGET_OBJECTS:kernel_${TARGET_CORE}>")
  endforeach()
endif ()

# add objects to the openblas lib
add_library(${OpenBLAS_LIBNAME} ${LA_SOURCES} ${LAPACKE_SOURCES} ${RELA_SOURCES} ${TARGET_OBJS} ${OpenBLAS_DEF_FILE})
target_include_directories(${OpenBLAS_LIBNAME} INTERFACE $<INSTALL_INTERFACE:include/openblas${SUFFIX64}>)

# Android needs to explicitly link against libm
if(ANDROID)
  target_link_libraries(${OpenBLAS_LIBNAME} m)
endif()

# Handle MSVC exports
if(MSVC AND BUILD_SHARED_LIBS)
  if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_LESS 3.4)
    include("${PROJECT_SOURCE_DIR}/cmake/export.cmake")
  else()
    # Creates verbose .def file (51KB vs 18KB)
    set_target_properties(${OpenBLAS_LIBNAME} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS true)
  endif()
endif()

# Set output for libopenblas
set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES LIBRARY_OUTPUT_NAME_DEBUG "${OpenBLAS_LIBNAME}_d")
set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES EXPORT_NAME "OpenBLAS")

foreach (OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
  string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )

  set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib/${OUTPUTCONFIG} )
  set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib/${OUTPUTCONFIG} )
  set_target_properties( ${OpenBLAS_LIBNAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR}/lib/${OUTPUTCONFIG} )
endforeach()

enable_testing()

if (USE_THREAD)
  # Add threading library to linker
  find_package(Threads)
  if (THREADS_HAVE_PTHREAD_ARG)
    set_property(TARGET ${OpenBLAS_LIBNAME} PROPERTY COMPILE_OPTIONS "-pthread")
    set_property(TARGET ${OpenBLAS_LIBNAME} PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
  endif()
  target_link_libraries(${OpenBLAS_LIBNAME} ${CMAKE_THREAD_LIBS_INIT})
endif()

#if (MSVC OR NOT NOFORTRAN)
if (NOT NO_CBLAS)
  # Broken without fortran on unix
  add_subdirectory(utest)
endif()

if (NOT MSVC AND NOT NOFORTRAN)
  # Build test and ctest
  add_subdirectory(test)
  if(NOT NO_CBLAS)
    add_subdirectory(ctest)
  endif()
  add_subdirectory(lapack-netlib/TESTING)
	  if (CPP_THREAD_SAFETY_TEST OR CPP_THREAD_SAFETY_GEMV)
	  	add_subdirectory(cpp_thread_test)
    	  endif()
endif()

set_target_properties(${OpenBLAS_LIBNAME} PROPERTIES
  VERSION ${OpenBLAS_MAJOR_VERSION}.${OpenBLAS_MINOR_VERSION}
  SOVERSION ${OpenBLAS_MAJOR_VERSION}
)

if (BUILD_SHARED_LIBS AND BUILD_RELAPACK)
  if (NOT MSVC)
    target_link_libraries(${OpenBLAS_LIBNAME} "-Wl,-allow-multiple-definition")
  else()
   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE")
  endif()
endif()

if (BUILD_SHARED_LIBS AND NOT ${SYMBOLPREFIX}${SYMBOLSUFFIX} STREQUAL "")
if (NOT DEFINED ARCH)
  set(ARCH_IN "x86_64")
else()
  set(ARCH_IN ${ARCH})
endif()

if (${CORE} STREQUAL "generic")
  set(ARCH_IN "GENERIC")
endif ()

if (NOT DEFINED EXPRECISION)
  set(EXPRECISION_IN 0)
else()
  set(EXPRECISION_IN ${EXPRECISION})
endif()

if (NOT DEFINED NO_CBLAS)
  set(NO_CBLAS_IN 0)
else()
  set(NO_CBLAS_IN ${NO_CBLAS})
endif()

if (NOT DEFINED NO_LAPACK)
  set(NO_LAPACK_IN 0)
else()
  set(NO_LAPACK_IN ${NO_LAPACK})
endif()

if (NOT DEFINED NO_LAPACKE)
  set(NO_LAPACKE_IN 0)
else()
  set(NO_LAPACKE_IN ${NO_LAPACKE})
endif()

if (NOT DEFINED NEED2UNDERSCORES)
  set(NEED2UNDERSCORES_IN 0)
else()
  set(NEED2UNDERSCORES_IN ${NEED2UNDERSCORES})
endif()

if (NOT DEFINED ONLY_CBLAS)
  set(ONLY_CBLAS_IN 0)
else()
  set(ONLY_CBLAS_IN ${ONLY_CBLAS})
endif()

if (NOT DEFINED BU)
  set(BU _)
endif()

if (NOT ${SYMBOLPREFIX} STREQUAL "")
message(STATUS "adding prefix ${SYMBOLPREFIX} to names of exported symbols in ${OpenBLAS_LIBNAME}")
endif()
if (NOT ${SYMBOLSUFFIX} STREQUAL "")
message(STATUS "adding suffix ${SYMBOLSUFFIX} to names of exported symbols in ${OpenBLAS_LIBNAME}")
endif()
	add_custom_command(TARGET ${OpenBLAS_LIBNAME} POST_BUILD
  	COMMAND perl  ${PROJECT_SOURCE_DIR}/exports/gensymbol "objcopy" "${ARCH}" "${BU}" "${EXPRECISION_IN}" "${NO_CBLAS_IN}" "${NO_LAPACK_IN}" "${NO_LAPACKE_IN}" "${NEED2UNDERSCORES_IN}" "${ONLY_CBLAS_IN}" \"${SYMBOLPREFIX}\" \"${SYMBOLSUFFIX}\" "${BUILD_LAPACK_DEPRECATED}" > ${PROJECT_BINARY_DIR}/objcopy.def
        COMMAND objcopy -v --redefine-syms ${PROJECT_BINARY_DIR}/objcopy.def  ${PROJECT_BINARY_DIR}/lib/lib${OpenBLAS_LIBNAME}.so
        COMMENT "renaming symbols"
        )
endif()


# Install project

# Install libraries
install(TARGETS ${OpenBLAS_LIBNAME}
	EXPORT "OpenBLAS${SUFFIX64}Targets"
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} )

# Install headers
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/openblas${SUFFIX64})
set(CMAKE_INSTALL_FULL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

message(STATUS "Generating openblas_config.h in ${CMAKE_INSTALL_INCLUDEDIR}")

set(OPENBLAS_CONFIG_H ${CMAKE_BINARY_DIR}/openblas_config.h)
file(WRITE  ${OPENBLAS_CONFIG_H} "#ifndef OPENBLAS_CONFIG_H\n")
file(APPEND ${OPENBLAS_CONFIG_H} "#define OPENBLAS_CONFIG_H\n")
file(STRINGS ${PROJECT_BINARY_DIR}/config.h __lines)
foreach(line ${__lines})
  string(REPLACE "#define " "" line ${line})
  file(APPEND ${OPENBLAS_CONFIG_H} "#define OPENBLAS_${line}\n")
endforeach()
file(APPEND ${OPENBLAS_CONFIG_H} "#define OPENBLAS_VERSION \"OpenBLAS ${OpenBLAS_VERSION}\"\n")
file(READ   ${CMAKE_CURRENT_SOURCE_DIR}/openblas_config_template.h OPENBLAS_CONFIG_TEMPLATE_H_CONTENTS)
file(APPEND ${OPENBLAS_CONFIG_H} "${OPENBLAS_CONFIG_TEMPLATE_H_CONTENTS}\n")
file(APPEND ${OPENBLAS_CONFIG_H} "#endif /* OPENBLAS_CONFIG_H */\n")
install (FILES ${OPENBLAS_CONFIG_H} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(NOT NOFORTRAN)
  message(STATUS "Generating f77blas.h in ${CMAKE_INSTALL_INCLUDEDIR}")

  set(F77BLAS_H ${CMAKE_BINARY_DIR}/generated/f77blas.h)
  file(WRITE  ${F77BLAS_H} "#ifndef OPENBLAS_F77BLAS_H\n")
  file(APPEND ${F77BLAS_H} "#define OPENBLAS_F77BLAS_H\n")
  file(APPEND ${F77BLAS_H} "#include \"openblas_config.h\"\n")
  file(READ ${CMAKE_CURRENT_SOURCE_DIR}/common_interface.h COMMON_INTERFACE_H_CONTENTS)
  file(APPEND ${F77BLAS_H} "${COMMON_INTERFACE_H_CONTENTS}\n")
  file(APPEND ${F77BLAS_H} "#endif")
  install (FILES ${F77BLAS_H} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

if(NOT NO_CBLAS)
	message (STATUS "Generating cblas.h in ${CMAKE_INSTALL_INCLUDEDIR}")
	set(CBLAS_H ${CMAKE_BINARY_DIR}/generated/cblas.h)
	file(READ ${CMAKE_CURRENT_SOURCE_DIR}/cblas.h CBLAS_H_CONTENTS)
	string(REPLACE "common" "openblas_config" CBLAS_H_CONTENTS_NEW "${CBLAS_H_CONTENTS}")
	if (NOT ${SYMBOLPREFIX} STREQUAL "")
	string(REPLACE " cblas" " ${SYMBOLPREFIX}cblas" CBLAS_H_CONTENTS	"${CBLAS_H_CONTENTS_NEW}")
	string(REPLACE " openblas" " ${SYMBOLPREFIX}openblas" CBLAS_H_CONTENTS_NEW	"${CBLAS_H_CONTENTS}")
	string (REPLACE " ${SYMBOLPREFIX}openblas_complex" " openblas_complex" CBLAS_H_CONTENTS	"${CBLAS_H_CONTENTS_NEW}")
	string(REPLACE " goto" " ${SYMBOLPREFIX}goto" CBLAS_H_CONTENTS_NEW "${CBLAS_H_CONTENTS}")
	endif()
	if (NOT ${SYMBOLSUFFIX} STREQUAL "")
	string(REGEX REPLACE "(cblas[^ (]*)" "\\1${SYMBOLSUFFIX}" CBLAS_H_CONTENTS	"${CBLAS_H_CONTENTS_NEW}")
	string(REGEX REPLACE "(openblas[^ (]*)" "\\1${SYMBOLSUFFIX}" CBLAS_H_CONTENTS_NEW "${CBLAS_H_CONTENTS}")
	string(REGEX REPLACE "(openblas_complex[^ ]*)${SYMBOLSUFFIX}" "\\1" CBLAS_H_CONTENTS	"${CBLAS_H_CONTENTS_NEW}")
	string(REGEX REPLACE "(goto[^ (]*)" "\\1${SYMBOLSUFFIX}" CBLAS_H_CONTENTS_NEW	"${CBLAS_H_CONTENTS}")
	endif()
	file(WRITE ${CBLAS_H} "${CBLAS_H_CONTENTS_NEW}")
	install (FILES ${CBLAS_H} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

if(NOT NO_LAPACKE)
	message (STATUS "Copying LAPACKE header files to ${CMAKE_INSTALL_INCLUDEDIR}")
	add_dependencies( ${OpenBLAS_LIBNAME} genlapacke)
	FILE(GLOB_RECURSE INCLUDE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lapack-netlib/LAPACKE/*.h")
	install (FILES ${INCLUDE_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

	ADD_CUSTOM_TARGET(genlapacke
	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/lapack-netlib/LAPACKE/include/lapacke_mangling_with_flags.h.in "${CMAKE_BINARY_DIR}/lapacke_mangling.h"
	)
	install (FILES ${CMAKE_BINARY_DIR}/lapacke_mangling.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/openblas${SUFFIX64})
endif()

# Install pkg-config files
configure_file(${PROJECT_SOURCE_DIR}/cmake/openblas.pc.in ${PROJECT_BINARY_DIR}/openblas${SUFFIX64}.pc @ONLY)
install (FILES ${PROJECT_BINARY_DIR}/openblas${SUFFIX64}.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/)


# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(PN OpenBLAS)
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}${SUFFIX64}")
configure_package_config_file(cmake/${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}${SUFFIX64}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                 VERSION ${${PN}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}${SUFFIX64}Config.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
        RENAME ${PN}${SUFFIX64}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}${SUFFIX64}Targets"
        NAMESPACE "${PN}${SUFFIX64}::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

