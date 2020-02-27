set(LAPACKE_INCLUDE lapacke.h lapack.h lapacke_config.h lapacke_utils.h)

file(COPY ${LAPACKE_INCLUDE} DESTINATION ${LAPACK_BINARY_DIR}/include)
