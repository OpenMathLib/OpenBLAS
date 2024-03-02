add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

lapack_major_version = 3 # soversion
lapack_minor_version = 12
lapack_patch_version = 0
lapack_version = f'@lapack_major_version@.@lapack_minor_version@.@lapack_patch_version@'

subdir('BLAS/SRC') # Defines blas and _allblas
subdir('INSTALL') # Defines _install_aux and _int_cpu_time
subdir('SRC') # Defines lapack
