##
## Author: Hank Anderson <hank@statease.com>
## Description: Ported from the OpenBLAS/c_check perl script.
##              This is triggered by prebuild.cmake and runs before any of the code is built.
##              Creates config.h and Makefile.conf.

# Convert CMake vars into the format that OpenBLAS expects
string(TOUPPER ${CMAKE_SYSTEM_NAME} HOST_OS)
if (${HOST_OS} STREQUAL "WINDOWS")
  set(HOST_OS WINNT)
endif ()

if(CMAKE_COMPILER_IS_GNUCC AND WIN32)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpmachine
              OUTPUT_VARIABLE OPENBLAS_GCC_TARGET_MACHINE
              OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(OPENBLAS_GCC_TARGET_MACHINE MATCHES "amd64|x86_64|AMD64")
      set(MINGW64 1)
    endif()
endif()

# Pretty thorough determination of arch. Add more if needed
if(CMAKE_CL_64 OR MINGW64)
  set(X86_64 1)
elseif(MINGW OR (MSVC AND NOT CMAKE_CROSSCOMPILING))
  set(X86 1)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc.*|power.*|Power.*")
  set(PPC 1)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "mips64.*")
  set(MIPS64 1)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "amd64.*|x86_64.*|AMD64.*")
  set(X86_64 1)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "i686.*|i386.*|x86.*|amd64.*|AMD64.*")
  set(X86 1)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*)")
  set(ARM 1)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*)")
  set(ARM64 1)
endif()

if (X86_64)
  set(ARCH "x86_64")
elseif(X86)
  set(ARCH "x86")
elseif(PPC)
  set(ARCH "power")
elseif(ARM)
  set(ARCH "arm")
elseif(ARM64)
  set(ARCH "arm64")
else()
  set(ARCH ${CMAKE_SYSTEM_PROCESSOR} CACHE STRING "Target Architecture")
endif ()

if (NOT BINARY)
  if (X86_64 OR ARM64 OR PPC OR MIPS64)
    set(BINARY 64)
  else ()
    set(BINARY 32)
  endif ()
endif()

if(BINARY EQUAL 64)
  set(BINARY64 1)
else()
  set(BINARY32 1)
endif()

