# TODO: Only Level 1 BLAS tests pass for now, the others need input files, see
# Makefile for more information
_blas_noinput_test_array = [#
  # ['Pretty name', 'binary_name', 'BlahTest.cpp']
]
_blas_input_test_array = [#
  # ['Pretty name', 'binary_name', 'BlahTest.cpp', 'inputfile.in']
]

if prec == 's' or build_single or build_all_prec
 _blas_noinput_test_array += [
   ['Test REAL Level 1 BLAS', 'xblat1s', 'sblat1.f'],
 ]
 # _blas_noinput_test_array += [
 #   ['Test REAL Level 2 BLAS', 'xblat2s', 'sblat2.f', 'sblat2.in'],
 #   ['Test REAL Level 3 BLAS', 'xblat3s', 'sblat3.f', 'sblat3.in'],
 # ]
endif

if prec == 'd' or build_double or build_all_prec
 _blas_noinput_test_array += [
   ['Test DOUBLE PRECISION Level 1 BLAS', 'xblat1d', 'dblat1.f'],
 ]
 # _blas_noinput_test_array += [
 #   ['Test DOUBLE PRECISION Level 2 BLAS', 'xblat2d', 'dblat2.f', 'dblat2.in'],
 #   ['Test DOUBLE PRECISION Level 3 BLAS', 'xblat3d', 'dblat3.f', 'dblat3.in'],
 # ]
endif


if build_complex or build_all_prec
 _blas_noinput_test_array += [
   ['Test COMPLEX Level 1 BLAS', 'xblat1c', 'cblat1.f'],
   # ['Test COMPLEX Level 2 BLAS', 'xblat2c', 'cblat2.f'],
   # ['Test COMPLEX Level 3 BLAS', 'xblat3c', 'cblat3.f'],
 ]
endif

if build_complex16 or build_all_prec
 _blas_noinput_test_array += [
   ['Test COMPLEX*16 Level 1 BLAS', 'xblat1z', 'zblat1.f'],
   # ['Test COMPLEX*16 Level 2 BLAS', 'xblat2z', 'zblat2.f'],
   # ['Test COMPLEX*16 Level 3 BLAS', 'xblat3z', 'zblat3.f'],
 ]
endif

foreach _test : _blas_noinput_test_array
  test(_test.get(0),
       executable(_test.get(1),
                  sources : [_test.get(2)],
                  link_with : blas,
                 ),
      )
endforeach

# foreach _test : _blas_input_test_array
#   test(_test.get(0),
#        executable(_test.get(1),
#                   sources : [_test.get(2)],
#                   link_with : blas,
#                  ),
#        args: ['<', _test.get(3)],
#        workdir : meson.source_root() + '/lapack-netlib/BLAS/TESTING/',
#       )
# endforeach
