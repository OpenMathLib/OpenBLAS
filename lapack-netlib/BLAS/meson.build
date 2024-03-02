_netlib_blas_name = get_option('netlib_blas_name')
_netlib_cblas_name = get_option('netlib_cblas_name')
_netlib_lapack_name = get_option('netlib_lapack_name')
_netlib_tmglib_name = get_option('netlib_tmglib_name')
_netlib_lapacke_name = get_option('netlib_lapacke_name')

subdir('SRC')

build_testing = get_option('build_testing')
if build_testing
  subdir('TESTING')
endif
