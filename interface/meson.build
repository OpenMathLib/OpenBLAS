_kinds = []
real_kinds = ['s', 'd']
if get_option('exprecision')
  real_kinds += 'q'
endif
complex_kinds = ['c', 'z', 'x']
extended_kinds = ['dx', 'bf16']
_rules = {
  's': ['FLOAT'],
  'd': ['DOUBLE'],
  'q': ['QUAD_PRECISION'],
  'c': ['COMPLEX'],
  'dx': ['XDOUBLE'],
  'bf19': ['BFLOAT16'],
}

blas_roots = []
blas1_roots = [
  # These don't exist as roots.
  # ismax amax isamax amin ismin isamin
  # TODO: Why is dsdot in sblas1objs!?
  # 'axpy', 'swap', 'copy', 'scal', 'dot', 'dsdot', 'asum', 'sum', 'nrm2', 'max',
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
