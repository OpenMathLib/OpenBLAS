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

sblas_extobs = []
_defs = []

blas_roots = []
blas1_roots = [
  # These don't exist as roots.
  # ismax amax isamax amin ismin isamin
  # TODO: Why is dsdot in sblas1objs!?
  'axpy', 'swap', 'copy', 'scal', 'dot', 'dsdot', 'asum', 'sum', 'nrm2', 'max',
  'rot', 'rotg', 'rotm', 'rotmg', 'axpby',
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
configurations = []
blas_roots += blas1_roots
blas_roots += blas2_roots
blas_roots += blas3_roots
_kinds += real_kinds
foreach blasi : blas_roots
  sfiles = blasi + '.c'
  # defs = []
  # if blasi in addl_srcs
  #   defs += addl_srcs[blasi]
  # endif
  foreach tkind : real_kinds
    name = tkind + blasi
    defs = []
    configurations += {'defs': defs, 'name': name, 'src': sfiles}
  endforeach
endforeach

# if get_option('bfloat16')
#   TODO: Handle bfloat16
#   sblas1_srcs += files('sbdot.c')
#   sblas2_srcs += files('sbgemv.c')
#   sblas3_srcs += files(['sbgemm.c', 'sbgemmt.c'])
# endif

_static_libs = []
foreach conf: configurations
    _static_libs += static_library(
        conf['name'],
        conf['src'],
        include_directories:  _inc,
        c_args: conf['defs'],
    )
endforeach

_interface = static_library('_interface',
                    link_whole: _static_libs)
