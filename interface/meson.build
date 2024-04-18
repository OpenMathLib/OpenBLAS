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
  # NOTE: q, qx, x, xq do not have cblas_ rules in the Makefile
  # NOTE: https://developer.arm.com/documentation/101004/2310/BLAS-Basic-Linear-Algebra-Subprograms/CBLAS-functions?lang=en
  # Level 1
  { 'base': '?asum', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'asum.c',
    'cblas': true,
  },
  { 'base': '?sum', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'sum.c',
    'cblas': true,
  },
  { 'base': '?amax', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'max.c', 'def': [ 'USE_ABS' ], 'undef': [ 'USE_MIN' ],
    'cblas': true,
  },
  { 'base': '?amin', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'max.c', 'def': [ 'USE_ABS', 'USE_MIN' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': '?max', '_types': ['s', 'd', 'q'],
    'fname': 'max.c', 'def': [ ], 'undef': [ 'USE_ABS', 'USE_MIN' ],
  },
  { 'base': '?min', '_types': ['s', 'd', 'q'],
    'fname': 'max.c', 'def': [ 'USE_MIN' ], 'undef': [ 'USE_ABS' ],
  },
  { 'base': '?axpy', '_types': ['s', 'd', 'q'],
    'fname': 'axpy.c',
    'cblas': true,
  },
  { 'base': '?axpy', '_types': ['c', 'z', 'x'],
    'fname': 'zaxpy.c',
    'cblas': true,
  },
  { 'base': '?axpyc', '_types': ['c', 'z', 'x'],
    'fname': 'zaxpy.c', 'def': [ 'CONJ' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': '?copy', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'copy.c',
    'cblas': true,
  },
  { 'base': '?dot', '_types': ['s', 'd', 'q'],
    'fname': 'dot.c',
    'cblas': true,
  },
  { 'base': '?dotu', '_types': ['c', 'z', 'x'],
    'fname': 'zdot.c', 'def': [ ], 'undef': [ 'CONJ' ],
    'cblas': true,
  },
  { 'base': '?dotc', '_types': ['c', 'z', 'x'],
    'fname': 'zdot.c', 'def': [ 'CONJ' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': 'cblas_?dotu_sub', '_types': ['c', 'z'],
    'fname': 'zdot.c', 'def': [ 'CBLAS', 'FORCE_USE_STACK' ], 'undef': [ 'CONJ' ],
    'cblas': false, # These don't non-cblas rules Makefile:1623:1627
  },
  { 'base': 'cblas_?dotc_sub', '_types': ['c', 'z'],
    'fname': 'zdot.c', 'def': [ 'CBLAS', 'FORCE_USE_STACK', 'CONJ' ], 'undef': [ ],
    'cblas': false, # These don't non-cblas rules Makefile:1623:1627
  },
  { 'base': '?nrm2', '_types': ['s', 'd', 'q', 'sc', 'dz', 'qx'],
    'fname': 'nrm2.c',
    'cblas': true, # q, qx don't have cblas_ rules Makefile:1635,1645
  },
  { 'base': '?rot', '_types': ['s', 'd', 'q'],
    'fname': 'rot.c',
    'cblas': true, # q doesn't have cblas_
  },
  { 'base': '?rot', '_types': ['cs', 'zd', 'xq'],
    'fname': 'zrot.c',
    'cblas': true,
  },
  { 'base': '?rotm', '_types': ['s', 'd', 'q'],
    'fname': 'rotm.c',
    'cblas': true,
  },
  { 'base': '?rotmg', '_types': ['s', 'd', 'q'],
    'fname': 'rotmg.c',
    'cblas': true,
  },
  { 'base': '?rotg', '_types': ['s', 'd', 'q'],
    'fname': 'rotg.c',
    'cblas': true,
  },
  { 'base': '?rotg', '_types': ['c', 'z', 'x'],
    'fname': 'zrotg.c',
    'cblas': true,
  },
  { 'base': '?scal', '_types': ['s', 'd', 'q'],
    'fname': 'scal.c',
    'cblas': true,
  },
  { 'base': '?scal', '_types': ['c', 'z', 'x'],
    'fname': 'zscal.c',
    'cblas': true, # x doesn't have cblas_
  },
  { 'base': '?scal', '_types': ['cs', 'zd', 'xq'],
    'fname': 'zscal.c', 'def': [ 'SSCAL' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': '?swap', '_types': ['s', 'd', 'q'],
    'fname': 'swap.c',
    'cblas': true,
  },
  { 'base': '?swap', '_types': ['c', 'z', 'x'],
    'fname': 'zswap.c',
    'cblas': true,
  },
  { 'base': '?dsdot', '_types': ['s'],
    'fname': 'sdsdot.c',
    'cblas': true,
  },
  { 'base': '?dsdot', '_types': [''],
    'fname': 'dsdot.c',
    'cblas': true,
  },
  # TODO: Handle BFLOAT16 Makefile 784:795, 1592:1603
  { 'base': 'i?max', '_types': ['s', 'd', 'q', 'c', 'z'], # NOTE: c,z only for cblas, Makefile:1526
    'fname': 'imax.c', 'def': [ ], 'undef': [ 'USE_ABS', 'USE_MIN' ],
    'cblas': true, # TODO: Maybe 'cblas_only': ['c', 'z']
  },
  { 'base': 'i?amax', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'imax.c', 'def': [ 'USE_ABS' ], 'undef': [ 'USE_MIN' ],
    'cblas': true,
  },
  { 'base': 'i?amin', '_types': ['s', 'd', 'c', 'z', 'x'],
    'fname': 'imax.c', 'def': [ 'USE_ABS', 'USE_MIN' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': 'i?min', '_types': ['s', 'd', 'q', 'c', 'z'], # NOTE: c,z only for cblas, Makefile:1532
    'fname': 'imax.c', 'def': [ 'USE_MIN' ], 'undef': [ 'USE_ABS' ],
    'cblas': true,
  },
  # Level 2
  { 'base': '?ger', '_types': ['s', 'd', 'q'],
    'fname': 'ger.c',
    'cblas': true,
  },
  { 'base': '?geru', '_types': ['c', 'z', 'x'],
    'fname': 'zger.c', 'def': [ ], 'undef': [ 'CONJ' ],
    'cblas': true,
  },
  { 'base': '?gerc', '_types': ['c', 'z', 'x'],
    'fname': 'zger.c', 'def': [ 'CONJ' ], 'undef': [ ],
    'cblas': true,
  },
  # TODO: Handle BFLOAT16 Makefile 941:944
  # TODO: Handle Netlib_gemv Makefile 946:958
  { 'base': '?gemv', '_types': ['s', 'd', 'q'],
    'fname': 'gemv.c',
    'cblas': true,
  },
  { 'base': '?gemv', '_types': ['c', 'z', 'x'],
    'fname': 'zgemv.c',
    'cblas': true,
  },
  { 'base': '?trmv', '_types': ['s', 'd', 'q'],
    'fname': 'trmv.c',
    'cblas': true,
  },
  { 'base': '?trmv', '_types': ['c', 'z', 'x'],
    'fname': 'ztrmv.c',
    'cblas': true,
  },
  { 'base': '?trsv', '_types': ['s', 'd', 'q'],
    'fname': 'trsv.c',
    'cblas': true,
  },
  { 'base': '?trsv', '_types': ['c', 'z', 'x'],
    'fname': 'ztrsv.c',
    'cblas': true,
  },
  # TODO: Handle NO_LAPACK
  { 'base': '?symv', '_types': ['s', 'd', 'q'],
    'fname': 'symv.c',
    'cblas': true,
  },
  { 'base': '?symv', '_types': ['c', 'z', 'x'],
    'fname': 'zsymv.c',
    'cblas': false, # NOTE: This oddly has no cblas set..
  },
  { 'base': '?syr', '_types': ['s', 'd', 'q'],
    'fname': 'syr.c',
    'cblas': true,
  },
  { 'base': '?syr', '_types': ['c', 'z', 'x'],
    'fname': 'zsyr.c',
    'cblas': false, # NOTE: This oddly has no cblas set..
  },
  { 'base': '?syr2', '_types': ['s', 'd', 'q'],
    'fname': 'syr2.c',
    'cblas': true,
  },
  { 'base': '?syr2', '_types': ['c', 'z', 'x'],
    'fname': 'zsyr2.c',
    'cblas': false, # NOTE: This oddly has no cblas set..
  },
  { 'base': '?gbmv', '_types': ['s', 'd', 'q'],
    'fname': 'gbmv.c',
    'cblas': true,
  },
  { 'base': '?gbmv', '_types': ['c', 'z', 'x'],
    'fname': 'zgbmv.c',
    'cblas': true,
  },
  { 'base': '?sbmv', '_types': ['s', 'd', 'q'],
    'fname': 'sbmv.c',
    'cblas': true,
  },
  { 'base': '?sbmv', '_types': ['c', 'z', 'x'],
    'fname': 'zsbmv.c',
    'cblas': false, # NOTE: This oddly has no cblas set..
  },
  { 'base': '?spmv', '_types': ['s', 'd', 'q'],
    'fname': 'spmv.c',
    'cblas': true,
  },
  # TODO: NO_LAPACK
  { 'base': '?spmv', '_types': ['c', 'z', 'x'],
    'fname': 'zspmv.c',
    'cblas': false, # NOTE: This oddly has no cblas set..
  },
  { 'base': '?spr', '_types': ['s', 'd', 'q'],
    'fname': 'spr.c',
    'cblas': true,
  },
  { 'base': '?spr', '_types': ['c', 'z', 'x'],
    'fname': 'zspr.c',
    'cblas': false, # NOTE: This oddly has no cblas set..
  },
  { 'base': '?spr2', '_types': ['s', 'd', 'q'],
    'fname': 'spr2.c',
    'cblas': true,
  },
  { 'base': '?spr2', '_types': ['c', 'z', 'x'],
    'fname': 'zspr2.c',
    'cblas': false, # NOTE: This oddly has no cblas set..
  },
  { 'base': '?tbmv', '_types': ['s', 'd', 'q'],
    'fname': 'tbmv.c',
    'cblas': true,
  },
  { 'base': '?tbmv', '_types': ['c', 'z', 'x'],
    'fname': 'ztbmv.c',
    'cblas': true,
  },
  { 'base': '?tbsv', '_types': ['s', 'd', 'q'],
    'fname': 'tbsv.c',
    'cblas': true,
  },
  { 'base': '?tbsv', '_types': ['c', 'z', 'x'],
    'fname': 'ztbsv.c',
    'cblas': true,
  },
  { 'base': '?tpsv', '_types': ['s', 'd', 'q'],
    'fname': 'tpsv.c',
    'cblas': true,
  },
  { 'base': '?tpsv', '_types': ['c', 'z', 'x'],
    'fname': 'ztpsv.c',
    'cblas': true,
  },
  { 'base': '?tpmv', '_types': ['s', 'd', 'q'],
    'fname': 'tpmv.c',
    'cblas': true,
  },
  { 'base': '?tpmv', '_types': ['c', 'z', 'x'],
    'fname': 'ztpmv.c',
    'cblas': true,
  },
  { 'base': '?hemv', '_types': ['c', 'z', 'x'],
    'fname': 'zhemv.c',
    'cblas': true,
  },
  { 'base': '?hbmv', '_types': ['c', 'z', 'x'],
    'fname': 'zhbmv.c',
    'cblas': true,
  },
  { 'base': '?her', '_types': ['c', 'z', 'x'],
    'fname': 'zher.c',
    'cblas': true,
  },
  { 'base': '?her2', '_types': ['c', 'z', 'x'],
    'fname': 'zher2.c',
    'cblas': true,
  },
  { 'base': '?hpmv', '_types': ['c', 'z', 'x'],
    'fname': 'zhpmv.c',
    'cblas': true,
  },
  { 'base': '?hpr', '_types': ['c', 'z', 'x'],
    'fname': 'zhpr.c',
    'cblas': true,
  },
  { 'base': '?hpr2', '_types': ['c', 'z', 'x'],
    'fname': 'zhpr2.c',
    'cblas': true,
  },
  # TODO: BUILD_BFLOAT16 here, Makefile:1303
  { 'base': '?gemm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': ['gemm.c', '../param.h'], # TODO: Will this work?
    'cblas': true,
  },
  { 'base': '?gemmt', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': ['gemmt.c', '../param.h'],
    'cblas': true,
  },
  { 'base': '?symm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'symm.c',
    'cblas': true,
  },
  { 'base': '?trmm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'trsm.c', 'def': [ 'TRMM' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': '?trsm', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'trsm.c',
    'cblas': true,
  },
  { 'base': '?syrk', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'syrk.c',
    'cblas': true,
  },
  { 'base': '?syr2k', '_types': ['s', 'd', 'q', 'c', 'z', 'x'],
    'fname': 'syr2k.c',
    'cblas': true,
  },
  { 'base': '?hemm', '_types': ['c', 'z', 'x'],
    'fname': 'symm.c', 'def': [ 'HEMM' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': '?herk', '_types': ['c', 'z', 'x'],
    'fname': 'syrk.c', 'def': [ 'HEMM' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': '?her2k', '_types': ['c', 'z', 'x'],
    'fname': 'syr2k.c', 'def': [ 'HEMM' ], 'undef': [ ],
    'cblas': true,
  },
  { 'base': '?gemm3m', '_types': ['c', 'z', 'x'],
    'fname': 'gemm.c', 'def': [ 'GEMM3M' ], 'undef': [ ],
    'cblas': true,
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
