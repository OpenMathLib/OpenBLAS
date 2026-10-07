message(STATUS "building experimental lapacke testsuite")
set (BUILD_DEFAULT_API ON)
# Precisions to test: those whose LAPACKE interface is built. A
# precision-pruned build tests only what it has.
set(LAPACKE_TEST_PRECISIONS)
foreach(_spec "s;SINGLE" "d;DOUBLE" "c;COMPLEX" "z;COMPLEX16")
  list(GET _spec 0 _prec)
  list(GET _spec 1 _build)
  if(BUILD_${_build})
    list(APPEND LAPACKE_TEST_PRECISIONS ${_prec})
  endif()
endforeach()
if(NOT LAPACKE_TEST_PRECISIONS)
message(STATUS "no lapacke tests selected")
  return()
endif()

add_subdirectory(wrappers)
add_subdirectory(high_level)

# The LIN tests, rerouted through LAPACKE.
#
# The test sources are rewritten so that calls to allowlisted routines
# resolve to the <NAME>_TEST wrappers above, which forward to LAPACKE; the
# unmodified numerical and error-exit checks then exercise LAPACKE. String
# replacements must stay off so that the SRNAMT names used by the error-exit
# tests keep matching the names the unrenamed LAPACK routines pass to the
# (interposed) testing XERBLA.
# include(ExtendedAPIHelpers)
include("${PROJECT_SOURCE_DIR}/cmake/ExtendedAPIHelpers.cmake")

# An empty SYMBOL_ALLOWLIST would rename every symbol, so pass the option
# only when there is something to route.
set(_allowlist_args)
if(LAPACKE_TEST_ALLOWLISTED_ROUTINES)
  list(APPEND _allowlist_args
    SYMBOL_ALLOWLIST ${LAPACKE_TEST_ALLOWLISTED_ROUTINES})
endif()

# Outputs live in the TESTING/lapacke subdirectory, where lapack_testing.py
# analyzes them as their own LAPACKE section. Only the work-level
# column-major flavor runs the error-exit tests: the row-major interface has
# different leading-dimension semantics (LDA >= N instead of LDA >= M), and
# the high-level interface ignores the caller's workspace so LWORK error
# exits cannot fire. The other flavors read a generated input with the
# error-exit tests disabled, and so does that one where the platform cannot
# run them at all.
set(_lapacke_test_out_dir ${PROJECT_BINARY_DIR}/TESTING/lapacke)
file(MAKE_DIRECTORY ${_lapacke_test_out_dir})

foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  string(TOUPPER "${_prec}" _prec_upper)
  if(NOT LAPACKE_TEST_LIN_${_prec_upper}_SOURCES)
    continue()
  endif()

  # Route the allowlisted calls to the wrappers. Both APIs start from this
  # same rewrite; the extended API then suffixes the result a second time.
  generate_suffixed_sources(xlintst${_prec}_lapacke
    LAPACKE_TEST_LIN_${_prec_upper}_SOURCES _lin_sources_test
    SUFFIX "_TEST"
    ${_allowlist_args}
    NO_STRING_REPLACEMENTS)
  lapack_nag_disable_constant_propagation(${_lin_sources_test})

  lapack_test_input(_test_input
	  "${PROJECT_SOURCE_DIR}/TESTING/${_prec}test.in"
    OUTPUT "${_lapacke_test_out_dir}/${_prec}test.in")
  lapack_test_input(_test_input_noerr
	  "${PROJECT_SOURCE_DIR}/TESTING/${_prec}test.in"
    DISABLE_ERROR_EXIT_TESTS
    OUTPUT "${_lapacke_test_out_dir}/${_prec}test_noerr.in")

  if(BUILD_DEFAULT_API)
    add_library(xlintst${_prec}_lapacke_core OBJECT ${_lin_sources_test})

    foreach(_layer work high)
      foreach(_layout cm rm)
        set(_driver xlintst${_prec}_${_layer}_${_layout})

        add_executable(${_driver}
          $<TARGET_OBJECTS:xlintst${_prec}_lapacke_core>)
        set_target_properties(${_driver} PROPERTIES LINKER_LANGUAGE Fortran)
        target_link_libraries(${_driver}
          lapacke_test_wrappers_${_layer}_${_layout}
          ${OpenBLAS_LIBNAME})

        set(_output "${_lapacke_test_out_dir}/${_prec}test_${_layer}_${_layout}.out")
        if(_layer STREQUAL "work" AND _layout STREQUAL "cm")
          set(_input "${_test_input}")
        else()
          set(_input "${_test_input_noerr}")
        endif()
        lapack_runtest_command(_test_command ${_driver}
          INPUT "${_input}" OUTPUT "${_output}")
        add_test(NAME LAPACKE-${_driver}_${_prec}test_in COMMAND ${_test_command})
      endforeach()
    endforeach()
  endif()

  if(BUILD_INDEX64_EXT_API)
    # Suffix every remaining symbol as well, so that the driver is ILP64
    # throughout and calls the _64 LAPACK. The routed calls become
    # <routine>_TEST_64, which is what the wrappers export under
    # LAPACKE_API64.
    generate_suffixed_sources(xlintst${_prec}_lapacke_TEST
      _lin_sources_test _lin_sources_test_64
      SUFFIX "_64")
    lapack_nag_disable_constant_propagation(${_lin_sources_test_64})

    add_library(xlintst${_prec}_lapacke_core_64 OBJECT ${_lin_sources_test_64})
    target_compile_options(xlintst${_prec}_lapacke_core_64
      PRIVATE ${FOPT_ILP64})

    foreach(_layer work high)
      foreach(_layout cm rm)
        set(_driver xlintst${_prec}_${_layer}_${_layout}_64)

        add_executable(${_driver}
          $<TARGET_OBJECTS:xlintst${_prec}_lapacke_core_64>)
        set_target_properties(${_driver} PROPERTIES LINKER_LANGUAGE Fortran)
        target_link_libraries(${_driver}
          lapacke_test_wrappers_${_layer}_${_layout}_64
	  ${TMGLIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${OpenBLAS_LIBNAME})


        set(_output "${_lapacke_test_out_dir}/${_prec}test_${_layer}_${_layout}_64.out")
        if(_layer STREQUAL "work" AND _layout STREQUAL "cm")
          set(_input "${_test_input}")
        else()
          set(_input "${_test_input_noerr}")
        endif()
        lapack_runtest_command(_test_command ${_driver}
          INPUT "${_input}" OUTPUT "${_output}")
        add_test(NAME LAPACKE-${_driver}_${_prec}test_in COMMAND ${_test_command})
      endforeach()
    endforeach()
  endif()
endforeach()
