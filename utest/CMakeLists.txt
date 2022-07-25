include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

if (MSVC AND "${CMAKE_C_COMPILER_ID}" MATCHES Clang)
  set(OpenBLAS_utest_src utest_main2.c)
else ()
  set(OpenBLAS_utest_src
    utest_main.c
    test_min.c
    test_amax.c
    test_ismin.c
    test_rotmg.c
    test_rot.c
    test_axpy.c
    test_dsdot.c
    test_dnrm2.c
    test_swap.c
  )
endif ()

# crashing on travis cl with an error code suggesting resource not found
if (NOT MSVC)
set(OpenBLAS_utest_src
  ${OpenBLAS_utest_src}
  test_dotu.c
  )
endif ()

# known to hang with the native Windows and Android threads
# FIXME needs checking if this works on any of the other platforms
if (OS_CYGWIN_NT OR OS_LINUX)
if (NOT USE_OPENMP)
set(OpenBLAS_utest_src
  ${OpenBLAS_utest_src}
  test_fork.c
  )
endif()
set(OpenBLAS_utest_src
  ${OpenBLAS_utest_src}
  test_post_fork.c
  )
endif()

if (NOT NO_LAPACK)
set(OpenBLAS_utest_src
  ${OpenBLAS_utest_src}
  test_potrs.c
  )
if (NOT NO_CBLAS AND NOT NO_LAPACKE)
set(OpenBLAS_utest_src
  ${OpenBLAS_utest_src}
  test_kernel_regress.c
  )
endif()
endif()

set(OpenBLAS_utest_bin openblas_utest)
add_executable(${OpenBLAS_utest_bin} ${OpenBLAS_utest_src})

target_link_libraries(${OpenBLAS_utest_bin} ${OpenBLAS_LIBNAME})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD" OR ${CMAKE_SYSTEM_NAME} MATCHES "QNX"  )
target_link_libraries(${OpenBLAS_utest_bin} m)
endif()

if (${CMAKE_SYSTEM_NAME} STREQUAL "WindowsStore")
set_target_properties( ${OpenBLAS_utest_bin} PROPERTIES COMPILE_DEFINITIONS "_CRT_SECURE_NO_WARNINGS")
endif()

#Set output for utest
set_target_properties( ${OpenBLAS_utest_bin} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
foreach (OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
  string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  set_target_properties( ${OpenBLAS_utest_bin} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

if (MSVC AND BUILD_SHARED_LIBS)
add_custom_command(TARGET ${OpenBLAS_utest_bin}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/lib/${CMAKE_CFG_INTDIR}/${OpenBLAS_LIBNAME}.dll ${CMAKE_CURRENT_BINARY_DIR}/.
)
endif()

add_test(${OpenBLAS_utest_bin} ${CMAKE_CURRENT_BINARY_DIR}/${OpenBLAS_utest_bin})
