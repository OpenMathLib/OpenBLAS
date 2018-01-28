include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

enable_language(Fortran)

set(OpenBLAS_Tests
  sblat1 sblat2 sblat3
  dblat1 dblat2 dblat3
  cblat1 cblat2 cblat3
  zblat1 zblat2 zblat3)

foreach(test_bin ${OpenBLAS_Tests})
add_executable(${test_bin} ${test_bin}.f)
target_link_libraries(${test_bin} ${OpenBLAS_LIBNAME})
endforeach()

set(float_types s d c z)
foreach(float_type ${float_types})
    string(TOUPPER ${float_type} float_type_upper)
    add_test(NAME "${float_type}blas1" 
      COMMAND "$<TARGET_FILE:${float_type}blat1>")
    add_test(NAME "${float_type}blas2" 
      COMMAND ${CMAKE_COMMAND}
        "-DTEST_COMMAND=$<TARGET_FILE:${float_type}blat2>"
        "-DTEST_INPUT=${PROJECT_SOURCE_DIR}/test/${float_type}blat2.dat"
        "-DTEST_FILE=${CMAKE_BINARY_DIR}/test/${float_type}blas2/${float_type_upper}BLAT2.SUMM"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/test_helper.cmake"
    )
    add_test(NAME "${float_type}blas3" 
      COMMAND ${CMAKE_COMMAND}
        "-DTEST_COMMAND=$<TARGET_FILE:${float_type}blat3>"
        "-DTEST_INPUT=${PROJECT_SOURCE_DIR}/test/${float_type}blat3.dat"
        "-DTEST_FILE=${CMAKE_BINARY_DIR}/test/${float_type}blas2/${float_type_upper}BLAT3.SUMM"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/test_helper.cmake"
    )
  endforeach()
