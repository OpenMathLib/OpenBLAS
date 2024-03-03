add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

lapack_major_version = 3 # soversion
lapack_minor_version = 12
lapack_patch_version = 0
lapack_version = f'@lapack_major_version@.@lapack_minor_version@.@lapack_patch_version@'

# Common variables
prec = get_option('realkind')
build_single = get_option('build_single')
build_double = get_option('build_double')
build_complex = get_option('build_complex')
build_complex16 = get_option('build_complex16')
build_all_prec = get_option('build_all_prec')
use_c_lapack = get_option('use_c_lapack')
use_xblas = get_option('use_xblas')

# Sub-directories
subdir('BLAS') # Defines netlib_blas and _allblas
subdir('INSTALL') # Defines _install_aux and _int_cpu_time
subdir('SRC') # Defines netlib_lapack
