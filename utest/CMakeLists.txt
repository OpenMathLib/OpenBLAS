include_directories(${PROJECT_SOURCE_DIR})

set(OpenBLAS_utest_src
  utest_main.c
  test_amax.c
  )

if (NOT NO_LAPACK)
set(OpenBLAS_utest_src
  ${OpenBLAS_utest_src}
  test_potrs.c
  )
endif()

set(OpenBLAS_utest_bin openblas_utest)
add_executable(${OpenBLAS_utest_bin} ${OpenBLAS_utest_src})

target_link_libraries(${OpenBLAS_utest_bin} ${OpenBLAS_LIBNAME})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
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

if (MSVC)
add_custom_command(TARGET ${OpenBLAS_utest_bin}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/lib/$<CONFIG>/${OpenBLAS_LIBNAME}.dll ${CMAKE_CURRENT_BINARY_DIR}/.
)
endif()

add_test(${OpenBLAS_utest_bin} ${CMAKE_CURRENT_BINARY_DIR}/${OpenBLAS_utest_bin})