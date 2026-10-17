# LAPACKE test wrapper libraries.
#
# The wrappers map <routine>_TEST Fortran symbols onto LAPACKE calls. The
# LAPACK testing programs are rewritten by the extended-API source generator
# (see TESTING/LIN/CMakeLists.txt) so that calls to allowlisted routines
# resolve to these wrappers. The same wrapper sources are compiled once per
# (layer, layout) combination:
#   layer:  work -- middle-level (_work) interface, caller workspace passes
#                   through (faithful WORK/LWORK semantics, error exits)
#           high -- high-level interface (NaN checks, internal allocation)
#   layout: cm   -- column-major LAPACKE interface (faithful passthrough)
#           rm   -- row-major interface through row-major shadow copies of
#                   all matrix arguments

set(LAPACKE_TEST_WRAPPER_SOURCES lapacke_test_wrappers.c)
foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  list(APPEND LAPACKE_TEST_WRAPPER_SOURCES lapacke_test_wrappers_${_prec}.c)
endforeach()

foreach(_layer work high)
  string(TOUPPER "${_layer}" _layer_upper)
  foreach(_layout cm rm)
    if(_layout STREQUAL "cm")
      set(_layout_define LAPACK_COL_MAJOR)
    else()
      set(_layout_define LAPACK_ROW_MAJOR)
    endif()

    if(BUILD_DEFAULT_API)
      set(_wrapper_lib lapacke_test_wrappers_${_layer}_${_layout})
      add_library(${_wrapper_lib} STATIC ${LAPACKE_TEST_WRAPPER_SOURCES})
      target_compile_definitions(${_wrapper_lib} PRIVATE
        LAPACKE_TEST_LAYOUT=${_layout_define}
        LAPACKE_TEST_LAYER=LAPACKE_TEST_LAYER_${_layer_upper})
      target_include_directories(${_wrapper_lib} PRIVATE
	      ${PROJECT_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/LAPACKE/include)
      target_link_libraries(${_wrapper_lib} PUBLIC ${LAPACKELIB})
    endif()

    if(BUILD_INDEX64_EXT_API)
      set(_wrapper_lib lapacke_test_wrappers_${_layer}_${_layout}_64)
      add_library(${_wrapper_lib} STATIC ${LAPACKE_TEST_WRAPPER_SOURCES})
      target_compile_definitions(${_wrapper_lib} PRIVATE
        LAPACKE_TEST_LAYOUT=${_layout_define}
        LAPACKE_TEST_LAYER=LAPACKE_TEST_LAYER_${_layer_upper}
        LAPACK_ILP64 LAPACKE_API64 WeirdNEC CBLAS_API64)
      target_include_directories(${_wrapper_lib} PRIVATE
        ${LAPACK_BINARY_DIR}/include ${LAPACK_SOURCE_DIR}/LAPACKE/include)
      target_link_libraries(${_wrapper_lib} PUBLIC ${LAPACKELIB})
    endif()
  endforeach()
endforeach()

# Generate the symbol allowlist: a routine is routed through LAPACKE if it
# has a LAPACKE implementation (a lapacke_<name>_work.c source) and a
# <NAME>_TEST wrapper in the LAPACKE test wrapper sources.
file(GLOB _lapacke_work_sources RELATIVE
	"${PROJECT_SOURCE_DIR}/lapack-netlib/LAPACKE/src"
	"${PROJECT_SOURCE_DIR}/lapack-netlib/LAPACKE/src/lapacke_*_work.c")

set(_lapacke_routines)
foreach(_work_source IN LISTS _lapacke_work_sources)
	message("trying ${_work_source} for work_sources")
  string(REGEX REPLACE "^lapacke_(.+)_work\\.c$" "\\1" _routine "${_work_source}")
  list(APPEND _lapacke_routines "${_routine}")
endforeach()
message (STATUS "list lapackeroutines is ${_lapacke_routines}")
set(_allowlisted_routines)
foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  set(_wrapper_source "${CMAKE_CURRENT_SOURCE_DIR}/lapacke_test_wrappers_${_prec}.c")
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_wrapper_source}")

  file(STRINGS "${_wrapper_source}" _wrapper_symbols REGEX "LAPACK_GLOBAL_SUFFIX\\([a-z0-9_]+_test,")
  foreach(_wrapper_symbol IN LISTS _wrapper_symbols)
    string(REGEX MATCH "LAPACK_GLOBAL_SUFFIX\\(([a-z0-9_]+)_test," _unused "${_wrapper_symbol}")
    set(_routine "${CMAKE_MATCH_1}")
    message(STATUS "matching routine ${_routine} against list")
    if(_routine IN_LIST _lapacke_routines)
    message(STATUS "adding wrapper for routine ${_routine}" )
      list(APPEND _allowlisted_routines "${_routine}")
    else()
      message(WARNING
        "LAPACKE test wrapper '${_routine}_test' has no matching LAPACKE"
        " routine and is not routed through LAPACKE")
    endif()
  endforeach()
endforeach()
list(REMOVE_DUPLICATES _allowlisted_routines)
list(SORT _allowlisted_routines)
list(LENGTH _allowlisted_routines _allowlisted_routine_count)
list(LENGTH _lapacke_routines _lapacke_routine_count)
message(STATUS
  "LAPACKE testing: ${_allowlisted_routine_count} of ${_lapacke_routine_count}"
  " routines routed through LAPACKE")

set(LAPACKE_TEST_ALLOWLISTED_ROUTINES ${_allowlisted_routines} PARENT_SCOPE)
