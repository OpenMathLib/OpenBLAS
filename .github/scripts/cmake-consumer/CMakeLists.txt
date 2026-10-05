cmake_minimum_required(VERSION 3.19)
project(openblas_cmake_consumer C)

find_package(OpenBLAS REQUIRED CONFIG)
# a repeated lookup must not try to redefine the imported target
find_package(OpenBLAS REQUIRED CONFIG)

get_target_property(implib OpenBLAS::OpenBLAS IMPORTED_IMPLIB)
get_filename_component(implib_name "${implib}" NAME)
if (NOT implib_name STREQUAL EXPECTED_IMPLIB)
  message(FATAL_ERROR "IMPORTED_IMPLIB is ${implib}, expected ${EXPECTED_IMPLIB}")
endif ()

add_executable(consumer main.c)
target_link_libraries(consumer PRIVATE OpenBLAS::OpenBLAS)
