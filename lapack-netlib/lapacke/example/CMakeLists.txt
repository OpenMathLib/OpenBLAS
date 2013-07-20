add_executable(xexample_DGESV_rowmajor example_DGESV_rowmajor.c)
add_executable(xexample_ZGESV_rowmajor example_ZGESV_rowmajor.c)

target_link_libraries(xexample_DGESV_rowmajor lapacke ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xexample_ZGESV_rowmajor lapacke ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_test(example_DGESV_rowmajor ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGESV_rowmajor)
add_test(example_ZGESV_rowmajor ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_ZGESV_rowmajor)
