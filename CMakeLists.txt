##
## Author: Hank Anderson <hank@statease.com>
## Copyright: (c) Stat-Ease, Inc.
## Created: 12/23/14
## Last Modified: 12/23/14
##

cmake_minimum_required(VERSION 2.8.4)
project(OpenBLAS)

# is this necessary? lapack-netlib has its own fortran checks in its CMakeLists.txt
#enable_language(Fortran)

message(WARNING "CMake support is experimental. This will not produce the same Makefiles that OpenBLAS ships with.")

include("${CMAKE_SOURCE_DIR}/cmake/system.cmake")


