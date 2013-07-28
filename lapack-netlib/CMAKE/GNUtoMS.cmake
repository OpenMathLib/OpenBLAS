
#=============================================================================
# GNUtoMS - CMake module for Windows import library conversion
# Copyright 2010-2011 Kitware, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of Kitware, Inc. nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# GNUtoMS works only for the GNU toolchain on Windows (MinGW and MSys).
set(GNUtoMS 0)
if(NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR NOT WIN32 OR CYGWIN)
  return()
endif()

# Locate auxiliary GNUtoMS files.
get_filename_component(GNUtoMS_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
set(GNUtoMS_DIR ${GNUtoMS_DIR}/GNUtoMS)

if(NOT CMAKE_SIZEOF_VOID_P)
  enable_language(C) # Find CMAKE_SIZEOF_VOID_P reliably.
endif()

# Find MS development environment setup script for this architecture.
if("${CMAKE_SIZEOF_VOID_P}" EQUAL 4)
  find_program(VCVARS32 NAMES vcvars32.bat
    PATHS
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\10.0\\Setup\\VC;ProductDir]/bin"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\9.0\\Setup\\VC;ProductDir]/bin"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0\\Setup\\VC;ProductDir]/bin"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\7.1\\Setup\\VC;ProductDir]/bin"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\6.0\\Setup\\Microsoft Visual C++;ProductDir]/bin"
    )
  set(GNUtoMS_ENV "${VCVARS32}")
  set(GNUtoMS_ARCH x86)
elseif("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
  find_program(VCVARSAMD64 NAMES vcvarsamd64.bat
    PATHS
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\10.0\\Setup\\VC;ProductDir]/bin/amd64"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\9.0\\Setup\\VC;ProductDir]/bin/amd64"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0\\Setup\\VC;ProductDir]/bin/amd64"
    )
  set(GNUtoMS_ENV "${VCVARSAMD64}")
  set(GNUtoMS_ARCH amd64)
endif()

if(GNUtoMS_ENV)
  set(GNUtoMS 1)

  # Create helper script to run lib.exe from MS environment.
  string(REPLACE "/" "\\" GNUtoMS_BAT "${GNUtoMS_ENV}")
  set(LIB ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/GNUtoMSlib.bat)
  configure_file(${GNUtoMS_DIR}/lib.bat.in ${LIB})

  # Teach CMake how to create a MS import library at link time.
  set(CMAKE_Fortran_CREATE_SHARED_LIBRARY
    "${CMAKE_Fortran_CREATE_SHARED_LIBRARY} -Wl,--output-def,<TARGET_NAME>.def"
    "<CMAKE_COMMAND> -Dlib=\"${LIB}\" -Ddef=\"<TARGET_NAME>.def\" -Ddll=\"<TARGET>\" -Dimp=\"<TARGET_IMPLIB>\" -P \"${GNUtoMS_DIR}/lib.cmake\""
    )
endif()
