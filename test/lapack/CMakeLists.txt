set(netlib_tests "${NETLIB_LAPACK_DIR}/TESTING")

# The Netlib drivers append the ILP64 suffix when constructing the target name.
if("${LIBNAMEPREFIX}openblas${LIBNAMESUFFIX}" STREQUAL "${OpenBLAS_LIBNAME}")
  set(SUFFIX64_UNDERSCORE "")
endif()
add_subdirectory("${netlib_tests}/LIN" LIN)
add_subdirectory("${netlib_tests}/EIG" EIG)

function(add_lapack_driver_test driver input output)
  set(test_name "LAPACK-${driver}_${input}_in")
  add_test(NAME "${test_name}" COMMAND "${CMAKE_COMMAND}"
    "-DTEST=$<TARGET_FILE:${driver}>"
    "-DINPUT=${netlib_tests}/${input}.in"
    "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/${output}.out"
    "-DINTDIR=."
    -P "${netlib_tests}/runtest.cmake")
  set_tests_properties("${test_name}" PROPERTIES FAIL_REGULAR_EXPRESSION "FATAL")
endfunction()

foreach(precision IN ITEMS s d c z)
  if(NOT TARGET xlintst${precision})
    continue()
  endif()
  add_lapack_driver_test(xlintst${precision} ${precision}test ${precision}test)
  add_lapack_driver_test(xlintstrf${precision} ${precision}test_rfp ${precision}test_rfp)
  foreach(case IN ITEMS nep sep se2 svd glm gqr gsv csd lse)
    add_lapack_driver_test(xeigtst${precision} ${case} ${precision}${case})
  endforeach()
  foreach(case IN ITEMS ec ed gg gd sb sg bal bak gbal gbak bb)
    add_lapack_driver_test(xeigtst${precision} ${precision}${case} ${precision}${case})
  endforeach()
  add_lapack_driver_test(xdmdeigtst${precision} ${precision}dmd ${precision}dmd)

  # The bundled EIG list uses the single-precision source for every GESVJ driver.
  set(driver test_${precision}gesvj_emptsw)
  if(NOT precision STREQUAL "s")
    set_property(TARGET ${driver} PROPERTY SOURCES "${netlib_tests}/EIG/${driver}.f90")
  endif()
  add_test(NAME LAPACK-${driver} COMMAND ${driver})
  add_test(NAME LAPACK-test_${precision}gges3_scaling COMMAND test_${precision}gges3_scaling)
endforeach()

if(TARGET xlintstds)
  add_lapack_driver_test(xlintstds dstest dstest)
endif()
if(TARGET xlintstzc)
  add_lapack_driver_test(xlintstzc zctest zctest)
endif()

if(PYTHONINTERP_FOUND)
  add_test(NAME LAPACK_Test_Summary
    COMMAND "${PYTHON_EXECUTABLE}" "${NETLIB_LAPACK_DIR}/lapack_testing.py"
      --dir "${CMAKE_CURRENT_BINARY_DIR}")
endif()
