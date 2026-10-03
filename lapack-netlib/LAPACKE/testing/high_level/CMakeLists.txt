# Dedicated tests for the high-level-only behavior the Fortran test suite
# cannot exercise: NaN checks actually detecting NaNs, and workspace /
# transposition allocation failures (via malloc failure injection).
#
# One test file per routine and precision in <x><name>.c, each defining one
# function via LAPACKE_TEST(<x><name>). There is one runner per precision,
# whose test table is generated from the list below; to add a test, add the
# file for every precision and list its routine name here.
set(LAPACKE_HIGH_LEVEL_ROUTINES
  gels geqrf getrf getri getrs lange matgen potrf potrs)

set(LAPACKE_HIGH_LEVEL_COMMON_SOURCES
  lapacke_test.c lapacke_test_alloc.c lapacke_test_nan_sweep.c)

foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  set(_sources ${LAPACKE_HIGH_LEVEL_COMMON_SOURCES} lapacke_test_matgen_${_prec}.c)

  # Test table for this precision, substituted into the runner main below.
  set(LAPACKE_HIGH_LEVEL_EXTERNS "")
  set(LAPACKE_HIGH_LEVEL_ENTRIES "")
  foreach(_routine IN LISTS LAPACKE_HIGH_LEVEL_ROUTINES)
    set(_test_name ${_prec}${_routine})
    list(APPEND _sources ${_test_name}.c)
    string(APPEND LAPACKE_HIGH_LEVEL_EXTERNS
      "extern void lapacke_test_fn_${_test_name}( void );\n")
    string(APPEND LAPACKE_HIGH_LEVEL_ENTRIES
      "    { \"${_test_name}\", lapacke_test_fn_${_test_name} },\n")
  endforeach()

  if(BUILD_DEFAULT_API)
    set(_runner ${_prec}lapacke_high_level)
    set(LAPACKE_HIGH_LEVEL_PROGRAM ${_runner})
    configure_file(lapacke_high_level_main.c.in
      ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c @ONLY)

    add_executable(${_runner} ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c ${_sources})
    target_include_directories(${_runner} PRIVATE
      ${LAPACK_BINARY_DIR}/include ${LAPACK_SOURCE_DIR}/LAPACKE/include)
    target_link_libraries(${_runner}
      ${OpenBLAS_LIBNAME})

    add_test(NAME LAPACKE-high-level-${_prec} COMMAND ${_runner})
  endif()

  if(BUILD_INDEX64_EXT_API)
    set(_runner ${_prec}lapacke_high_level_64)
    set(LAPACKE_HIGH_LEVEL_PROGRAM ${_runner})
    configure_file(lapacke_high_level_main.c.in
      ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c @ONLY)

    add_executable(${_runner} ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c ${_sources})
    target_compile_definitions(${_runner} PRIVATE
      LAPACK_ILP64 LAPACKE_API64 WeirdNEC CBLAS_API64)
    target_include_directories(${_runner} PRIVATE
      ${LAPACK_BINARY_DIR}/include ${LAPACK_SOURCE_DIR}/LAPACKE/include)
    target_link_libraries(${_runner}
      ${LAPACKELIB} ${TMGLIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

    add_test(NAME LAPACKE-high-level-${_prec}_64 COMMAND ${_runner})
  endif()
endforeach()
