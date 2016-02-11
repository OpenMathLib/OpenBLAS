include_directories(${CMAKE_SOURCE_DIR})

set(OpenBLAS_utest_src
  utest_main.c
  test_amax.c
  )

set(OpenBLAS_utest_bin openblas_utest)
add_executable(${OpenBLAS_utest_bin} ${OpenBLAS_utest_src})
target_link_libraries(${OpenBLAS_utest_bin} ${OpenBLAS_LIBNAME})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(${OpenBLAS_utest_bin} m)
endif()

add_test(${OpenBLAS_utest_bin} ${CMAKE_CURRENT_BINARY_DIR}/${OpenBLAS_utest_bin})