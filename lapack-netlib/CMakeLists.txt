cmake_minimum_required(VERSION 3.2)

project(LAPACK Fortran C)

set(LAPACK_MAJOR_VERSION 3)
set(LAPACK_MINOR_VERSION 11)
set(LAPACK_PATCH_VERSION 0)
set(
  LAPACK_VERSION
  ${LAPACK_MAJOR_VERSION}.${LAPACK_MINOR_VERSION}.${LAPACK_PATCH_VERSION}
  )

# Allow setting a prefix for the library names
set(CMAKE_STATIC_LIBRARY_PREFIX "lib${LIBRARY_PREFIX}")
set(CMAKE_SHARED_LIBRARY_PREFIX "lib${LIBRARY_PREFIX}")

# Add the CMake directory for custom CMake modules
set(CMAKE_MODULE_PATH "${LAPACK_SOURCE_DIR}/CMAKE" ${CMAKE_MODULE_PATH})

# Export all symbols on Windows when building shared libraries
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "Coverage")
endif()

# Coverage
set(_is_coverage_build 0)
set(_msg "Checking if build type is 'Coverage'")
message(STATUS "${_msg}")
if(NOT CMAKE_CONFIGURATION_TYPES)
  string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type_lc)
  if(${_build_type_lc} STREQUAL "coverage")
    set(_is_coverage_build 1)
  endif()
endif()
message(STATUS "${_msg}: ${_is_coverage_build}")

if(_is_coverage_build)
  message(STATUS "Adding coverage")
  find_package(codecov)
endif()

# By default test Fortran compiler complex abs and complex division
option(TEST_FORTRAN_COMPILER "Test Fortran compiler complex abs and complex division" OFF)
if( TEST_FORTRAN_COMPILER )

  add_executable( test_zcomplexabs ${LAPACK_SOURCE_DIR}/INSTALL/test_zcomplexabs.f )
  add_custom_target( run_test_zcomplexabs
    COMMAND test_zcomplexabs 2> test_zcomplexabs.err
    WORKING_DIRECTORY ${LAPACK_BINARY_DIR}/INSTALL
    COMMENT "Running test_zcomplexabs in ${LAPACK_BINARY_DIR}/INSTALL with stderr: test_zcomplexabs.err"
    SOURCES ${LAPACK_SOURCE_DIR}/INSTALL/test_zcomplexabs.f )

  add_executable( test_zcomplexdiv ${LAPACK_SOURCE_DIR}/INSTALL/test_zcomplexdiv.f )
  add_custom_target( run_test_zcomplexdiv
    COMMAND test_zcomplexdiv 2> test_zcomplexdiv.err
    WORKING_DIRECTORY ${LAPACK_BINARY_DIR}/INSTALL
    COMMENT "Running test_zcomplexdiv in ${LAPACK_BINARY_DIR}/INSTALL with stderr: test_zcomplexdiv.err"
    SOURCES ${LAPACK_SOURCE_DIR}/INSTALL/test_zcomplexdiv.f )

  add_executable( test_zcomplexmult ${LAPACK_SOURCE_DIR}/INSTALL/test_zcomplexmult.f )
  add_custom_target( run_test_zcomplexmult
    COMMAND test_zcomplexmult 2> test_zcomplexmult.err
    WORKING_DIRECTORY ${LAPACK_BINARY_DIR}/INSTALL
    COMMENT "Running test_zcomplexmult in ${LAPACK_BINARY_DIR}/INSTALL with stderr: test_zcomplexmult.err"
    SOURCES ${LAPACK_SOURCE_DIR}/INSTALL/test_zcomplexmult.f )

  add_executable( test_zminMax ${LAPACK_SOURCE_DIR}/INSTALL/test_zminMax.f )
  add_custom_target( run_test_zminMax
    COMMAND test_zminMax 2> test_zminMax.err
    WORKING_DIRECTORY ${LAPACK_BINARY_DIR}/INSTALL
    COMMENT "Running test_zminMax in ${LAPACK_BINARY_DIR}/INSTALL with stderr: test_zminMax.err"
    SOURCES ${LAPACK_SOURCE_DIR}/INSTALL/test_zminMax.f )
    
endif()

# By default static library
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# By default build index32 library
option(BUILD_INDEX64 "Build Index-64 API libraries" OFF)
if(BUILD_INDEX64)
  set(BLASLIB "blas64")
  set(CBLASLIB "cblas64")
  set(LAPACKLIB "lapack64")
  set(LAPACKELIB "lapacke64")
  set(TMGLIB "tmglib64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWeirdNEC -DLAPACK_ILP64 -DHAVE_LAPACK_CONFIG_H")
  set(FORTRAN_ILP TRUE)
else()
  set(BLASLIB "blas")
  set(CBLASLIB "cblas")
  set(LAPACKLIB "lapack")
  set(LAPACKELIB "lapacke")
  set(TMGLIB "tmglib")
endif()

include(GNUInstallDirs)

# Updated OSX RPATH settings
# In response to CMake 3.0 generating warnings regarding policy CMP0042,
# the OSX RPATH settings have been updated per recommendations found
# in the CMake Wiki:
#  http://www.cmake.org/Wiki/CMake_RPATH_handling#Mac_OS_X_and_the_RPATH
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${CMAKE_INSTALL_FULL_LIBDIR} isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR})
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()


# Configure the warning and code coverage suppression file
configure_file(
  "${LAPACK_SOURCE_DIR}/CTestCustom.cmake.in"
  "${LAPACK_BINARY_DIR}/CTestCustom.cmake"
  @ONLY
)

include(PreventInSourceBuilds)
include(PreventInBuildInstalls)

# Check if recursive flag exists
include(CheckFortranCompilerFlag)
if(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
  check_fortran_compiler_flag("-Mrecursive" _MrecursiveFlag)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  check_fortran_compiler_flag("-frecursive" _frecursiveFlag)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  check_fortran_compiler_flag("-recursive" _recursiveFlag)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  check_fortran_compiler_flag("-qrecur" _qrecurFlag)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  check_fortran_compiler_flag("-recursive" _recursiveFlag)
else()
  message(WARNING "Fortran local arrays should be allocated on the stack."
    " Please use a compiler which guarantees that feature."
    " See https://github.com/Reference-LAPACK/lapack/pull/188 and references therein.")
endif()

# Add recursive flag
if(_MrecursiveFlag)
  string(REGEX MATCH "-Mrecursive" output_test <string> "${CMAKE_Fortran_FLAGS}")
  if(NOT output_test)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mrecursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()
elseif(_frecursiveFlag)
  string(REGEX MATCH "-frecursive" output_test <string> "${CMAKE_Fortran_FLAGS}")
  if(NOT output_test)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()
elseif(_recursiveFlag)
  string(REGEX MATCH "-recursive" output_test <string> "${CMAKE_Fortran_FLAGS}")
  if(NOT output_test)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()
elseif(_qrecurFlag)
  string(REGEX MATCH "-qrecur" output_test <string> "${CMAKE_Fortran_FLAGS}")
  if(NOT output_test)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qrecur"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()
endif()

if(UNIX)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model strict")
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnosave -qstrict")
  endif()
# Delete libmtsk in linking sequence for Sun/Oracle Fortran Compiler.
# This library is not present in the Sun package SolarisStudio12.3-linux-x86-bin
  string(REPLACE \;mtsk\; \; CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL Compaq)
  if(WIN32)
    if(CMAKE_GENERATOR STREQUAL "NMake Makefiles")
      get_filename_component(CMAKE_Fortran_COMPILER_CMDNAM ${CMAKE_Fortran_COMPILER} NAME_WE)
      message(STATUS "Using Compaq Fortran compiler with command name ${CMAKE_Fortran_COMPILER_CMDNAM}")
      set(cmd ${CMAKE_Fortran_COMPILER_CMDNAM})
      string(TOLOWER "${cmd}" cmdlc)
      if(cmdlc STREQUAL "df")
        message(STATUS "Assume the Compaq Visual Fortran Compiler is being used")
        set(CMAKE_Fortran_USE_RESPONSE_FILE_FOR_OBJECTS 1)
        set(CMAKE_Fortran_USE_RESPONSE_FILE_FOR_INCLUDES 1)
        #This is a workaround that is needed to avoid forward-slashes in the
        #filenames listed in response files from incorrectly being interpreted as
        #introducing compiler command options
        if(${BUILD_SHARED_LIBS})
          message(FATAL_ERROR "Making of shared libraries with CVF has not been tested.")
        endif()
        set(str "NMake version 9 or later should be used. NMake version 6.0 which is\n")
        set(str "${str}   included with the CVF distribution fails to build Lapack because\n")
        set(str "${str}   the number of source files exceeds the limit for NMake v6.0\n")
        message(STATUS ${str})
        set(CMAKE_Fortran_LINK_EXECUTABLE "LINK /out:<TARGET> <LINK_FLAGS> <LINK_LIBRARIES> <OBJECTS>")
      endif()
    endif()
  endif()
endif()

# Add option to enable flat namespace for symbol resolution on macOS
if(APPLE)
  option(USE_FLAT_NAMESPACE "Use flat namespaces for symbol resolution during build and runtime." OFF)

  if(USE_FLAT_NAMESPACE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-flat_namespace")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,-flat_namespace")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-flat_namespace")
  else()
    if(BUILD_SHARED_LIBS AND BUILD_TESTING)
      message(WARNING
        "LAPACK test suite might fail with shared libraries and the default two-level namespace. "
        "Disable shared libraries or enable flat namespace for symbol resolution via -DUSE_FLAT_NAMESPACE=ON.")
    endif()
  endif()
endif()

# --------------------------------------------------
set(LAPACK_INSTALL_EXPORT_NAME ${LAPACKLIB}-targets)

macro(lapack_install_library lib)
  install(TARGETS ${lib}
    EXPORT ${LAPACK_INSTALL_EXPORT_NAME}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
  )
endmacro()

set(PKG_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# --------------------------------------------------
# Testing
option(BUILD_TESTING "Build tests" ${_is_coverage_build})
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")

if(BUILD_TESTING)
  set(_msg "Looking for Python3 needed for summary tests")
  message(STATUS "${_msg}")
  # find_package(PythonInterp 3) cannot be used because /usr/bin/python may be
  # a Python2 interpreter.
  find_program(PYTHON_EXECUTABLE python3)
  if(PYTHON_EXECUTABLE)
    message(STATUS "${_msg} - found")
  else()
    message(STATUS "${_msg} - not found (skipping summary tests)")
  endif()
endif()

# --------------------------------------------------
# Organize output files.  On Windows this also keeps .dll files next
# to the .exe files that need them, making tests easy to run.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${LAPACK_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LAPACK_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LAPACK_BINARY_DIR}/lib)

# --------------------------------------------------
# Check for any necessary platform specific compiler flags
include(CheckLAPACKCompilerFlags)
CheckLAPACKCompilerFlags()

# --------------------------------------------------
# Check second function

include(CheckTimeFunction)
set(TIME_FUNC NONE)
CHECK_TIME_FUNCTION(NONE TIME_FUNC)
CHECK_TIME_FUNCTION(INT_CPU_TIME TIME_FUNC)
CHECK_TIME_FUNCTION(EXT_ETIME TIME_FUNC)
CHECK_TIME_FUNCTION(EXT_ETIME_ TIME_FUNC)
CHECK_TIME_FUNCTION(INT_ETIME TIME_FUNC)
message(STATUS "--> Will use second_${TIME_FUNC}.f and dsecnd_${TIME_FUNC}.f as timing function.")

set(SECOND_SRC ${LAPACK_SOURCE_DIR}/INSTALL/second_${TIME_FUNC}.f)
set(DSECOND_SRC ${LAPACK_SOURCE_DIR}/INSTALL/dsecnd_${TIME_FUNC}.f)

# deprecated LAPACK and LAPACKE routines
option(BUILD_DEPRECATED "Build deprecated routines" OFF)
message(STATUS "Build deprecated routines: ${BUILD_DEPRECATED}")

# --------------------------------------------------
# Precision to build
# By default all precisions are generated
option(BUILD_SINGLE "Build single precision real" ON)
option(BUILD_DOUBLE "Build double precision real" ON)
option(BUILD_COMPLEX "Build single precision complex" ON)
option(BUILD_COMPLEX16 "Build double precision complex" ON)
message(STATUS "Build single precision real: ${BUILD_SINGLE}")
message(STATUS "Build double precision real: ${BUILD_DOUBLE}")
message(STATUS "Build single precision complex: ${BUILD_COMPLEX}")
message(STATUS "Build double precision complex: ${BUILD_COMPLEX16}")

if(NOT (BUILD_SINGLE OR BUILD_DOUBLE OR BUILD_COMPLEX OR BUILD_COMPLEX16))
  message(FATAL_ERROR "Nothing to build, no precision selected.
  Please enable at least one of these:
  BUILD_SINGLE, BUILD_DOUBLE, BUILD_COMPLEX, BUILD_COMPLEX16.")
endif()


# --------------------------------------------------
# Subdirectories that need to be processed
option(USE_OPTIMIZED_BLAS "Whether or not to use an optimized BLAS library instead of included netlib BLAS" OFF)

# Check the usage of the user provided BLAS libraries
if(BLAS_LIBRARIES)
  include(CheckFortranFunctionExists)
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
  CHECK_FORTRAN_FUNCTION_EXISTS("dgemm" BLAS_FOUND)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(BLAS_FOUND)
    message(STATUS "--> BLAS supplied by user is WORKING, will use ${BLAS_LIBRARIES}.")
  else()
    message(ERROR "--> BLAS supplied by user is not WORKING, CANNOT USE ${BLAS_LIBRARIES}.")
    message(ERROR "-->     Will use REFERENCE BLAS (by default)")
    message(ERROR "-->     Or Correct your BLAS_LIBRARIES entry ")
    message(ERROR "-->     Or Consider checking USE_OPTIMIZED_BLAS")
  endif()

# User did not provide a BLAS Library but specified to search for one
elseif(USE_OPTIMIZED_BLAS)
  find_package(BLAS)
endif()

# Neither user specified or optimized BLAS libraries can be used
if(NOT BLAS_FOUND)
  message(STATUS "Using supplied NETLIB BLAS implementation")
  add_subdirectory(BLAS)
  set(BLAS_LIBRARIES ${BLASLIB})
else()
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    CACHE STRING "Linker flags for executables" FORCE)
  set(CMAKE_MODULE_LINKER_FLAGS
    "${CMAKE_MODULE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    CACHE STRING "Linker flags for modules" FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    CACHE STRING "Linker flags for shared libs" FORCE)
endif()


# --------------------------------------------------
# CBLAS
option(CBLAS "Build CBLAS" OFF)

if(CBLAS)
  add_subdirectory(CBLAS)
endif()

# --------------------------------------------------
# XBLAS

option(USE_XBLAS "Build extended precision (needs XBLAS)" OFF)
if(USE_XBLAS)
  find_library(XBLAS_LIBRARY NAMES xblas)
endif()

option(USE_OPTIMIZED_LAPACK "Whether or not to use an optimized LAPACK library instead of included netlib LAPACK" OFF)

# --------------------------------------------------
# LAPACK
# User did not provide a LAPACK Library but specified to search for one
if(USE_OPTIMIZED_LAPACK)
  find_package(LAPACK)
endif()

# Check the usage of the user provided or automatically found LAPACK libraries
if(LAPACK_LIBRARIES)
  include(CheckFortranFunctionExists)
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
  # Check if new routine of 3.4.0 is in LAPACK_LIBRARIES
  CHECK_FORTRAN_FUNCTION_EXISTS("dgeqrt" LATESTLAPACK_FOUND)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(LATESTLAPACK_FOUND)
    message(STATUS "--> LAPACK supplied by user is WORKING, will use ${LAPACK_LIBRARIES}.")
  else()
    message(ERROR "--> LAPACK supplied by user is not WORKING or is older than LAPACK 3.4.0, CANNOT USE ${LAPACK_LIBRARIES}.")
    message(ERROR "-->     Will use REFERENCE LAPACK (by default)")
    message(ERROR "-->     Or Correct your LAPACK_LIBRARIES entry ")
    message(ERROR "-->     Or Consider checking USE_OPTIMIZED_LAPACK")
  endif()
endif()

# Neither user specified or optimized LAPACK libraries can be used
if(NOT LATESTLAPACK_FOUND)
  message(STATUS "Using supplied NETLIB LAPACK implementation")
  set(LAPACK_LIBRARIES ${LAPACKLIB})
  add_subdirectory(SRC)
else()
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}"
    CACHE STRING "Linker flags for executables" FORCE)
  set(CMAKE_MODULE_LINKER_FLAGS
    "${CMAKE_MODULE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}"
    CACHE STRING "Linker flags for modules" FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}"
    CACHE STRING "Linker flags for shared libs" FORCE)
endif()

if(BUILD_TESTING)
  add_subdirectory(TESTING)
endif()

# --------------------------------------------------
# LAPACKE
option(LAPACKE "Build LAPACKE" OFF)

# LAPACKE has also the interface to some routines from tmglib,
# if LAPACKE_WITH_TMG is selected, we need to add those routines to LAPACKE
option(LAPACKE_WITH_TMG "Build LAPACKE with tmglib routines" OFF)
if(LAPACKE_WITH_TMG)
  set(LAPACKE ON)
endif()

# TMGLIB
# Cache export target
set(LAPACK_INSTALL_EXPORT_NAME_CACHE ${LAPACK_INSTALL_EXPORT_NAME})
if(BUILD_TESTING OR LAPACKE_WITH_TMG)
  if(LATESTLAPACK_FOUND AND LAPACKE_WITH_TMG)
    set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
    # Check if dlatms (part of tmg) is found
    CHECK_FORTRAN_FUNCTION_EXISTS("dlatms" LAPACK_WITH_TMGLIB_FOUND)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(NOT LAPACK_WITH_TMGLIB_FOUND)
      # Build and install TMG as part of LAPACKE targets (as opposed to LAPACK
      # targets)
      set(LAPACK_INSTALL_EXPORT_NAME ${LAPACKELIB}-targets)
    endif()
  endif()
  add_subdirectory(TESTING/MATGEN)
endif()
# Reset export target
set(LAPACK_INSTALL_EXPORT_NAME ${LAPACK_INSTALL_EXPORT_NAME_CACHE})
unset(LAPACK_INSTALL_EXPORT_NAME_CACHE)

if(LAPACKE)
  add_subdirectory(LAPACKE)
endif()


#-------------------------------------
# BLAS++ / LAPACK++
option(BLAS++ "Build BLAS++" OFF)
option(LAPACK++ "Build LAPACK++" OFF)


function(_display_cpp_implementation_msg name)
  string(TOLOWER ${name} name_lc)
  message(STATUS "${name}++ enable")
  message(STATUS "----------------")
  message(STATUS "Thank you for your interest in ${name}++, a newly developed C++ API for ${name} library")
  message(STATUS "The objective of ${name}++ is to provide a convenient, performance oriented API for development in the C++ language, that, for the most part, preserves established conventions, while, at the same time, takes advantages of modern C++ features, such as: namespaces, templates, exceptions, etc.")
  message(STATUS "For support ${name}++ related question, please email: slate-user@icl.utk.edu")
  message(STATUS "----------------")
endfunction()
if (BLAS++)
  _display_cpp_implementation_msg("BLAS")
  include(ExternalProject)
  ExternalProject_Add(blaspp
	  URL https://bitbucket.org/icl/blaspp/downloads/blaspp-2020.10.02.tar.gz
	  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}/lib LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} ${PROJECT_BINARY_DIR}/blaspp-prefix/src/blaspp
	  BUILD_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib LIB_SUFFIX="" make
	  INSTALL_COMMAND make PREFIX=${PROJECT_BINARY_DIR} LIB_SUFFIX="" install
	  )
  ExternalProject_Add_StepDependencies(blaspp build ${BLAS_LIBRARIES})
endif()
if (LAPACK++)
	message (STATUS "linking lapack++ against ${LAPACK_LIBRARIES}")
  _display_cpp_implementation_msg("LAPACK")
  include(ExternalProject)
  if (BUILD_SHARED_LIBS)
  ExternalProject_Add(lapackpp
	  URL https://bitbucket.org/icl/lapackpp/downloads/lapackpp-2020.10.02.tar.gz
	  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}/lib LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}  ${PROJECT_BINARY_DIR}/lapackpp-prefix/src/lapackpp
	  BUILD_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib LIB_SUFFIX="" make
	  INSTALL_COMMAND make PREFIX=${PROJECT_BINARY_DIR} LIB_SUFFIX="" install
	  )
  else ()
# FIXME this does not really work as the libraries list gets converted to a semicolon-separated list somewhere in the lapack++ build files
  ExternalProject_Add(lapackpp
	  URL https://bitbucket.org/icl/lapackpp/downloads/lapackpp-2020.10.02.tar.gz
	  CONFIGURE_COMMAND env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}/lib LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DLAPACK_LIBRARIES="${PROJECT_BINARY_DIR}/lib/liblapack.a -lgfortran"  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}  ${PROJECT_BINARY_DIR}/lapackpp-prefix/src/lapackpp
	  BUILD_COMMAND env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib LIB_SUFFIX="" make
	  INSTALL_COMMAND make PREFIX=${PROJECT_BINARY_DIR} LIB_SUFFIX="" install
  )
  endif()
  ExternalProject_Add_StepDependencies(lapackpp build blaspp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# --------------------------------------------------
# CPACK Packaging

set(CPACK_PACKAGE_NAME "LAPACK")
set(CPACK_PACKAGE_VENDOR "University of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LAPACK- Linear Algebra Package")
set(CPACK_PACKAGE_VERSION_MAJOR ${LAPACK_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${LAPACK_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${LAPACK_PATCH_VERSION})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_MONOLITHIC_INSTALL ON)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "LAPACK")
if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backslashes.
  set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\http://icl.cs.utk.edu/lapack-forum")
  set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.netlib.org/lapack")
  set(CPACK_NSIS_CONTACT "lapack@eecs.utk.edu")
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_NSIS_DISPLAY_NAME "LAPACK-${LAPACK_VERSION}")
  set(CPACK_PACKAGE_RELOCATABLE "true")
else()
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR TGZ)
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "lapack-${LAPACK_VERSION}")
  set(CPACK_SOURCE_IGNORE_FILES ~$ .svn ${CPACK_SOURCE_IGNORE_FILES})
endif()
include(CPack)


# --------------------------------------------------

if(NOT BLAS_FOUND)
  set(ALL_TARGETS ${ALL_TARGETS} ${BLASLIB})
endif()

if(NOT LATESTLAPACK_FOUND)
  set(ALL_TARGETS ${ALL_TARGETS} ${LAPACKLIB})
endif()

# Export lapack targets, not including lapacke, from the
# install tree, if any.
set(_lapack_config_install_guard_target "")
if(ALL_TARGETS)
  install(EXPORT ${LAPACKLIB}-targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LAPACKLIB}-${LAPACK_VERSION}
    COMPONENT Development
    )

  # Choose one of the lapack targets to use as a guard for
  # lapack-config.cmake to load targets from the install tree.
  list(GET ALL_TARGETS 0 _lapack_config_install_guard_target)
endif()

# Include cblas in targets exported from the build tree.
if(CBLAS)
  set(ALL_TARGETS ${ALL_TARGETS} ${CBLASLIB})
endif()

# Include lapacke in targets exported from the build tree.
if(LAPACKE)
  set(ALL_TARGETS ${ALL_TARGETS} ${LAPACKELIB})
endif()

if(NOT LAPACK_WITH_TMGLIB_FOUND AND LAPACKE_WITH_TMG)
  set(ALL_TARGETS ${ALL_TARGETS} ${TMGLIB})
endif()

# Export lapack and lapacke targets from the build tree, if any.
set(_lapack_config_build_guard_target "")
if(ALL_TARGETS)
  export(TARGETS ${ALL_TARGETS} FILE ${LAPACKLIB}-targets.cmake)

  # Choose one of the lapack or lapacke targets to use as a guard
  # for lapack-config.cmake to load targets from the build tree.
  list(GET ALL_TARGETS 0 _lapack_config_build_guard_target)
endif()

configure_file(${LAPACK_SOURCE_DIR}/CMAKE/lapack-config-build.cmake.in
  ${LAPACK_BINARY_DIR}/${LAPACKLIB}-config.cmake @ONLY)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lapack.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKLIB}.pc @ONLY)
  install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKLIB}.pc
  DESTINATION ${PKG_CONFIG_DIR}
  COMPONENT Development
  )

configure_file(${LAPACK_SOURCE_DIR}/CMAKE/lapack-config-install.cmake.in
  ${LAPACK_BINARY_DIR}/CMakeFiles/${LAPACKLIB}-config.cmake @ONLY)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${LAPACK_BINARY_DIR}/${LAPACKLIB}-config-version.cmake
  VERSION ${LAPACK_VERSION}
  COMPATIBILITY SameMajorVersion
  )

install(FILES
  ${LAPACK_BINARY_DIR}/CMakeFiles/${LAPACKLIB}-config.cmake
  ${LAPACK_BINARY_DIR}/${LAPACKLIB}-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LAPACKLIB}-${LAPACK_VERSION}
  COMPONENT Development
  )
if (LAPACK++)
  install(
  DIRECTORY "${LAPACK_BINARY_DIR}/lib/"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  FILES_MATCHING REGEX "liblapackpp.(a|so)$"
  )
  install(
    DIRECTORY "${PROJECT_BINARY_DIR}/lapackpp-prefix/src/lapackpp/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING REGEX "\\.(h|hh)$"
  )
  write_basic_package_version_file(
    "lapackppConfigVersion.cmake"
    VERSION 2020.10.02
    COMPATIBILITY AnyNewerVersion
  )
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/lib/lapackpp/lapackppConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/lib/lapackpp/lapackppConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/"
  )

endif()
if (BLAS++)
  write_basic_package_version_file(
    "blasppConfigVersion.cmake"
    VERSION 2020.10.02
    COMPATIBILITY AnyNewerVersion
  )
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/lib/blaspp/blasppConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/lib/blaspp/blasppConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/"
  )
  install(
  DIRECTORY "${LAPACK_BINARY_DIR}/lib/"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  FILES_MATCHING REGEX "libblaspp.(a|so)$"
  )
  install(
    DIRECTORY "${PROJECT_BINARY_DIR}/blaspp-prefix/src/blaspp/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING REGEX "\\.(h|hh)$"
  )
endif()

# --------------------------------------------------
# Generate MAN and/or HTML Documentation
option(BUILD_HTML_DOCUMENTATION "Create and install the HTML based API
documentation (requires Doxygen) - command: make html" OFF)
option(BUILD_MAN_DOCUMENTATION "Create and install the MAN based documentation (requires Doxygen) - command: make man" OFF)
message(STATUS "Build html documentation: ${BUILD_HTML_DOCUMENTATION}")
message(STATUS "Build man documentation: ${BUILD_MAN_DOCUMENTATION}")

if(BUILD_HTML_DOCUMENTATION OR BUILD_MAN_DOCUMENTATION)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
         message(WARNING "Doxygen is needed to build the documentation.")

    else()

    set(DOXYGEN_PROJECT_BRIEF "LAPACK: Linear Algebra PACKage")
    set(DOXYGEN_PROJECT_NUMBER ${LAPACK_VERSION})
    set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/DOCS)
    set(PROJECT_LOGO ${CMAKE_CURRENT_SOURCE_DIR}/DOCS/lapack.png)
    set(DOXYGEN_OPTIMIZE_FOR_FORTRAN YES)
    set(DOXYGEN_SOURCE_BROWSER YES)
    set(DISTRIBUTE_GROUP_DOC YES)
    set(DOXYGEN_CREATE_SUBDIRS YES)
    set(DOXYGEN_SEPARATE_MEMBER_PAGES YES)
    set(DOXYGEN_EXTRACT_ALL YES)
    set(DOXYGEN_FILE_PATTERNS "*.f;*.c;*.h")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_INTERACTIVE_SVG YES)
    set(DOXYGEN_QUIET YES)
    set(DOXYGEN_WARNINGS NO)
    set(DOXYGEN_GENERATE_HTML NO)
    set(DOXYGEN_GENERATE_MAN NO)


    if (BUILD_HTML_DOCUMENTATION)
        set(DOXYGEN_GENERATE_HTML YES)
        set(DOXYGEN_HTML_OUTPUT explore-html)
        set(DOXYGEN_INLINE_SOURCES YES)
        set(DOXYGEN_CALL_GRAPH YES)
        set(DOXYGEN_CALLER_GRAPH YES)

        doxygen_add_docs(
            html
            ${PROJECT_SOURCE_DIR}
            COMMENT "Generating html LAPACK documentation (it will take some time... time to grab a coffee)"
        )
    endif()
    if (BUILD_MAN_DOCUMENTATION)
        set(DOXYGEN_GENERATE_MAN YES)
        set(DOXYGEN_EXCLUDE SRC/VARIANTS)
        set(DOXYGEN_MAN_LINKS YES)
        set(DOXYGEN_INLINE_SOURCES NO)
        set(DOXYGEN_CALL_GRAPH NO)
        set(DOXYGEN_CALLER_GRAPH NO)

        doxygen_add_docs(
            man
            ${PROJECT_SOURCE_DIR}
            COMMENT "Generating man LAPACK documentation"
        )
    endif()

    endif()
endif()
