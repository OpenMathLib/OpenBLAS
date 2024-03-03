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
 _blas_input_test_array += [
   ['Test REAL Level 2 BLAS', 'xblat2s', 'sblat2.f', 'sblat2.in'],
   ['Test REAL Level 3 BLAS', 'xblat3s', 'sblat3.f', 'sblat3.in'],
 ]
endif

if prec == 'd' or build_double or build_all_prec
 _blas_noinput_test_array += [
   ['Test DOUBLE PRECISION Level 1 BLAS', 'xblat1d', 'dblat1.f'],
 ]
 _blas_input_test_array += [
   ['Test DOUBLE PRECISION Level 2 BLAS', 'xblat2d', 'dblat2.f', 'dblat2.in'],
   ['Test DOUBLE PRECISION Level 3 BLAS', 'xblat3d', 'dblat3.f', 'dblat3.in'],
 ]
endif


if prec == 'c' or build_complex or build_all_prec
 _blas_noinput_test_array += [
   ['Test COMPLEX Level 1 BLAS', 'xblat1c', 'cblat1.f'],
 ]
 _blas_input_test_array += [
   ['Test COMPLEX Level 2 BLAS', 'xblat2c', 'cblat2.f', 'cblat2.in'],
   ['Test COMPLEX Level 3 BLAS', 'xblat3c', 'cblat3.f', 'cblat3.in'],
 ]
endif

if prec == 'z' or build_complex16 or build_all_prec
 _blas_noinput_test_array += [
   ['Test COMPLEX*16 Level 1 BLAS', 'xblat1z', 'zblat1.f'],
 ]
 _blas_input_test_array += [
   ['Test COMPLEX*16 Level 2 BLAS', 'xblat2z', 'zblat2.f', 'zblat2.in'],
   ['Test COMPLEX*16 Level 3 BLAS', 'xblat3z', 'zblat3.f', 'zblat3.in'],
 ]
endif

foreach _test : _blas_noinput_test_array
  test(_test.get(0),
       executable(_test.get(1),
                  sources : _test.get(2),
                  link_with : netlib_blas,
                 ),
      )
endforeach

fortran_test_runner = executable('run_fortran_test',
                                 sources: ['run_fortran.c'],
                                 install: false)

# NOTE: For the tests to pass the executables need to be compiled first
foreach _test : _blas_input_test_array
  executable(_test.get(1),
             sources : _test.get(2),
             link_with : netlib_blas,
            )
  test_exe = meson.current_build_dir() / _test.get(1)
  input_file = meson.source_root() + '/lapack-netlib/BLAS/TESTING/' + _test.get(3)
  test(_test.get(0), fortran_test_runner,
       args: [test_exe, input_file],
       workdir: meson.current_build_dir())
endforeach
