_kinds = []
real_kinds = ['s', 'd']
if get_option('exprecision')
  real_kinds += 'q'
endif
complex_kinds = ['c', 'z', 'x']
extended_kinds = ['dx', 'bf16']

blas_roots = []
blas1_roots = [
  # These don't exist as roots.
  # ismax amax isamax amin ismin isamin
  # TODO: Why is dsdot in sblas1objs!?
  'axpy', 'swap', 'copy', 'scal', 'dot', 'dsdot', 'asum', 'sum', 'nrm2', 'max',
  'rot'#, 'rotg', 'rotm', 'rotmg', 'axpby',
]
blas2_roots = [
  'gemv', 'ger', 'trsv', 'trmv', 'symv', 'syr', 'syr2', 'gbmv', 'sbmv',
  'spmv', 'spr', 'spr2', 'tbsv', 'tbmv', 'tpsv', 'tpmv',
]
blas3_roots = [
  # TODO: trmm is trsm with a flag
  'gemm', 'symm', 'trsm', 'syrk', 'syr2k', 'omatcopy', 'imatcopy', 'geadd',
  'gemmt'
]

_blas_roots = [
  # Level 1
  { 'base': '?asum', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'asum.c',
  },
  { 'base': '?sum', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'sum.c',
  },
  { 'base': '?amax', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'max.c', 'def': [ 'USE_ABS' ], 'undef': [ 'USE_MIN' ],
  },
  { 'base': '?amin', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'max.c', 'def': [ 'USE_ABS', 'USE_MIN' ], 'undef': [ ],
  },
  { 'base': '?max', '_types': ['s', 'd', 'q'],
    'fname': 'max.c', 'def': [ ], 'undef': [ 'USE_ABS', 'USE_MIN' ],
  },
  { 'base': '?min', '_types': ['s', 'd', 'q'],
    'fname': 'max.c', 'def': [ 'USE_MIN' ], 'undef': [ 'USE_ABS' ],
  },
  { 'base': '?axpy', '_types': ['s', 'd', 'q'],
    'fname': 'axpy.c',
  },
  { 'base': '?axpy', '_types': ['c', 'z', 'x'],
    'fname': 'zaxpy.c',
  },
  { 'base': '?axpyc', '_types': ['c', 'z', 'x'],
    'fname': 'zaxpy.c', 'def': [ 'CONJ' ], 'undef': [ ],
  },
  { 'base': '?copy', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'copy.c',
  },
  { 'base': '?dot', '_types': ['s', 'd', 'q'],
    'fname': 'dot.c',
  },
  { 'base': '?dotu', '_types': ['c', 'z', 'x'],
    'fname': 'zdot.c', 'def': [ ], 'undef': [ 'CONJ' ],
  },
  { 'base': '?dotc', '_types': ['c', 'z', 'x'],
    'fname': 'zdot.c', 'def': [ 'CONJ' ], 'undef': [ ],
  },
  { 'base': '?nrm2', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'nrm2.c',
  },
  { 'base': '?rot', '_types': ['s', 'd', 'q'],
    'fname': 'rot.c',
  },
  { 'base': '?rot', '_types': ['cs', 'zd', 'xq'],
    'fname': 'zrot.c',
  },
  { 'base': '?rotm', '_types': ['s', 'd', 'q'],
    'fname': 'rotm.c',
  },
  { 'base': '?rotmg', '_types': ['s', 'd', 'q'],
    'fname': 'rotmg.c',
  },
  { 'base': '?rotg', '_types': ['s', 'd', 'q'],
    'fname': 'rotg.c',
  },
  { 'base': '?rotg', '_types': ['c', 'z', 'x'],
    'fname': 'zrotg.c',
  },
  { 'base': '?scal', '_types': ['s', 'd', 'q'],
    'fname': 'scal.c',
  },
  { 'base': '?scal', '_types': ['c', 'z', 'x'],
    'fname': 'zscal.c',
  },
  { 'base': '?scal', '_types': ['cs', 'zd', 'xq'],
    'fname': 'zscal.c', 'def': [ 'SSCAL' ], 'undef': [ ],
  },
  { 'base': '?swap', '_types': ['s', 'd', 'q'],
    'fname': 'swap.c',
  },
  { 'base': '?swap', '_types': ['c', 'z', 'x'],
    'fname': 'zswap.c',
  },
  { 'base': '?dsdot', '_types': ['s'],
    'fname': 'sdsdot.c',
  },
  { 'base': '?dsdot', '_types': [''],
    'fname': 'dsdot.c',
  },
  # TODO: Handle BFLOAT16 Makefile 784:795
  { 'base': 'i?max', '_types': ['s', 'd', 'q'],
    'fname': 'imax.c', 'def': [ ], 'undef': [ 'USE_ABS', 'USE_MIN' ],
  },
  { 'base': 'i?amax', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'imax.c', 'def': [ 'USE_ABS' ], 'undef': [ 'USE_MIN' ],
  },
  { 'base': 'i?amin', '_types': ['s', 'd', 'c', 'z', 'x'],
    'fname': 'imax.c', 'def': [ 'USE_ABS', 'USE_MIN' ], 'undef': [ ],
  },
  { 'base': 'i?min', '_types': ['s', 'd', 'q'],
    'fname': 'imax.c', 'def': [ 'USE_MIN' ], 'undef': [ 'USE_ABS' ],
  },
  # Level 2
  { 'base': '?ger', '_types': ['s', 'd', 'q'],
    'fname': 'ger.c',
  },
  { 'base': '?geru', '_types': ['c', 'z', 'x'],
    'fname': 'zger.c', 'def': [ ], 'undef': [ 'CONJ' ],
  },
  { 'base': '?gerc', '_types': ['c', 'z', 'x'],
    'fname': 'zger.c', 'def': [ 'CONJ' ], 'undef': [ ],
  },
  # TODO: Handle BFLOAT16 Makefile 941:944
  # TODO: Handle Netlib_gemv Makefile 946:958
  { 'base': '?gemv', '_types': ['s', 'd', 'q'],
    'fname': 'gemv.c',
  },
  { 'base': '?gemv', '_types': ['c', 'z', 'x'],
    'fname': 'zgemv.c',
  },
  { 'base': '?trmv', '_types': ['s', 'd', 'q'],
    'fname': 'trmv.c',
  },
  { 'base': '?trmv', '_types': ['c', 'z', 'x'],
    'fname': 'ztrmv.c',
  },
  { 'base': '?trsv', '_types': ['s', 'd', 'q'],
    'fname': 'trsv.c',
  },
  { 'base': '?trsv', '_types': ['c', 'z', 'x'],
    'fname': 'ztrsv.c',
  },
  # TODO: Handle NO_LAPACK
  { 'base': '?symv', '_types': ['s', 'd', 'q'],
    'fname': 'symv.c',
  },
  { 'base': '?symv', '_types': ['c', 'z', 'x'],
    'fname': 'zsymv.c',
  },
  { 'base': '?syr', '_types': ['s', 'd', 'q'],
    'fname': 'syr.c',
  },
  { 'base': '?syr', '_types': ['c', 'z', 'x'],
    'fname': 'zsyr.c',
  },
  { 'base': '?syr2', '_types': ['s', 'd', 'q'],
    'fname': 'syr2.c',
  },
  { 'base': '?syr2', '_types': ['c', 'z', 'x'],
    'fname': 'zsyr2.c',
  },
  { 'base': '?gbmv', '_types': ['s', 'd', 'q'],
    'fname': 'gbmv.c',
  },
  { 'base': '?gbmv', '_types': ['c', 'z', 'x'],
    'fname': 'zgbmv.c',
  },
  { 'base': '?sbmv', '_types': ['s', 'd', 'q'],
    'fname': 'sbmv.c',
  },
  { 'base': '?sbmv', '_types': ['c', 'z', 'x'],
    'fname': 'zsbmv.c',
  },
  { 'base': '?spmv', '_types': ['s', 'd', 'q'],
    'fname': 'spmv.c',
  },
  # TODO: NO_LAPACK
  { 'base': '?spmv', '_types': ['c', 'z', 'x'],
    'fname': 'zspmv.c',
  },
  { 'base': '?spr', '_types': ['s', 'd', 'q'],
    'fname': 'spr.c',
  },
  { 'base': '?spr', '_types': ['c', 'z', 'x'],
    'fname': 'zspr.c',
  },
  { 'base': '?spr2', '_types': ['s', 'd', 'q'],
    'fname': 'spr2.c',
  },
  { 'base': '?spr2', '_types': ['c', 'z', 'x'],
    'fname': 'zspr2.c',
  },
  { 'base': '?tbmv', '_types': ['s', 'd', 'q'],
    'fname': 'tbmv.c',
  },
  { 'base': '?tbmv', '_types': ['c', 'z', 'x'],
    'fname': 'ztbmv.c',
  },
  { 'base': '?tbsv', '_types': ['s', 'd', 'q'],
    'fname': 'tbsv.c',
  },
  { 'base': '?tbsv', '_types': ['c', 'z', 'x'],
    'fname': 'ztbsv.c',
  },
  { 'base': '?tpsv', '_types': ['s', 'd', 'q'],
    'fname': 'tpsv.c',
  },
  { 'base': '?tpsv', '_types': ['c', 'z', 'x'],
    'fname': 'ztpsv.c',
  },
  { 'base': '?tpmv', '_types': ['s', 'd', 'q'],
    'fname': 'tpmv.c',
  },
  { 'base': '?tpmv', '_types': ['c', 'z', 'x'],
    'fname': 'ztpmv.c',
  },
  { 'base': '?hemv', '_types': ['c', 'z', 'x'],
    'fname': 'zhemv.c',
  },
  { 'base': '?hbmv', '_types': ['c', 'z', 'x'],
    'fname': 'zhbmv.c',
  },
  { 'base': '?her', '_types': ['c', 'z', 'x'],
    'fname': 'zher.c',
  },
  { 'base': '?her2', '_types': ['c', 'z', 'x'],
    'fname': 'zher2.c',
  },
  { 'base': '?hpmv', '_types': ['c', 'z', 'x'],
    'fname': 'zhpmv.c',
  },
  { 'base': '?hpr', '_types': ['c', 'z', 'x'],
    'fname': 'zhpr.c',
  },
  { 'base': '?hpr2', '_types': ['c', 'z', 'x'],
    'fname': 'zhpr2.c',
  },
  # TODO: BUILD_BFLOAT16 here, Makefile:1303
  { 'base': '?gemm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': ['gemm.c', '../param.h'], # TODO: Will this work?
  },
  { 'base': '?gemmt', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': ['gemmt.c', '../param.h'],
  },
  { 'base': '?symm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'symm.c',
  },
  { 'base': '?trmm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'trsm.c', 'def': [ 'TRMM' ], 'undef': [ ],
  },
  { 'base': '?trsm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'trsm.c',
  },
  { 'base': '?syrk', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'syrk.c',
  },
  { 'base': '?syr2k', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'syr2k.c',
  },
  { 'base': '?hemm', '_types': ['c', 'z', 'x'],
    'fname': 'symm.c', 'def': [ 'HEMM' ], 'undef': [ ],
  },
  { 'base': '?herk', '_types': ['c', 'z', 'x'],
    'fname': 'syrk.c', 'def': [ 'HEMM' ], 'undef': [ ],
  },
  { 'base': '?her2k', '_types': ['c', 'z', 'x'],
    'fname': 'syr2k.c', 'def': [ 'HEMM' ], 'undef': [ ],
  },
  { 'base': '?gemm3m', '_types': ['c', 'z', 'x'],
    'fname': 'gemm.c', 'def': [ 'GEMM3M' ], 'undef': [ ],
  },
  { 'base': '?symm3m', '_types': ['c', 'z', 'x'],
    'fname': 'symm.c', 'def': [ 'GEMM3M' ], 'undef': [ ],
  },
  { 'base': '?hemm3m', '_types': ['c', 'z', 'x'],
    'fname': 'symm.c', 'def': [ 'HEMM', 'GEMM3M' ], 'undef': [ ],
  },
]

# addl_srcs = {
#   'symm': '-DHEMM'
# }

# Generated
# NOTE: Remember to check sourcesets as a simpler mechanism
# TODO: Generated the imax and max via use_abs use_min
interface_confs = []
blas_roots += blas1_roots
blas_roots += blas2_roots
blas_roots += blas3_roots
_kinds += real_kinds
foreach root : blas1_roots
  fname = root + '.c'
  defs = []
  # if root in addl_srcs
  #   defs += addl_srcs[root]
  # endif
  foreach prec : real_kinds
    _defs = []
    name = prec + fname
    if prec == 'd'
      _defs += ['-UCOMPLEX', '-DDOUBLE']
    elif prec == 's'
      _defs += ['-UCOMPLEX', '-UDOUBLE']
    elif prec == 'cs'
      _defs += ['-DCOMPLEX', '-UDOUBLE']
    endif
    interface_confs += {'defs': _defs,
                        'name': prec + root,
                        'src': fname}
    interface_confs += {'defs': _defs + '-DCBLAS',
                        'name': 'cblas_' + prec + root,
                        'src': fname}
  endforeach
endforeach
message(interface_confs)

# if get_option('bfloat16')
#   TODO: Handle bfloat16
#   sblas1_srcs += files('sbdot.c')
#   sblas2_srcs += files('sbgemv.c')
#   sblas3_srcs += files(['sbgemm.c', 'sbgemmt.c'])
# endif

_static_libs = []
foreach conf: interface_confs
  sym_name = conf['name']
  sym_underscored = f'@sym_name@_'
    _static_libs += static_library(
        sym_name,
        conf['src'],
        include_directories:  _inc,
        c_args: [
          conf['defs'],
          _cargs,
          f'-DASMNAME=@sym_name@',
          f'-DASMFNAME=@sym_underscored@',
          f'-DNAME=@sym_underscored@',
          f'-DCNAME=@sym_name@',
          f'-DCHAR_NAME="@sym_underscored@"',
          f'-DCHAR_CNAME="@sym_name@"',
        ],
    )
endforeach

_interface = static_library('_interface',
                            c_args: _cargs,
                            link_whole: _static_libs)
