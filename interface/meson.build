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
interface_confs = []
blas_roots += blas1_roots
blas_roots += blas2_roots
blas_roots += blas3_roots
_kinds += real_kinds
foreach root : blas_roots
  fname = root + '.c'
  defs = []
  # if root in addl_srcs
  #   defs += addl_srcs[root]
  # endif
  foreach prec : real_kinds
    name = prec + fname
    interface_confs += {'defs': defs, 'name': name, 'src': fname}
  endforeach
endforeach

# if get_option('bfloat16')
#   TODO: Handle bfloat16
#   sblas1_srcs += files('sbdot.c')
#   sblas2_srcs += files('sbgemv.c')
#   sblas3_srcs += files(['sbgemm.c', 'sbgemmt.c'])
# endif

interface_confs += {
  'name': 'cblas_drot',
  'defs': ['-DCBLAS'],
  'src': files('rot.c')
}

_static_libs = []
# foreach conf: interface_confs
#     _static_libs += static_library(
#         conf['name'],
#         conf['src'],
#         include_directories:  _inc,
#         c_args: conf['defs'],
#     )
# endforeach
_static_libs += static_library(
    'cblas_drot',
    'rot.c',
    include_directories:  _inc,
    c_args: [_cargs,
        '-DSMP_SERVER',
        '-DNO_WARMUP',
        '-DMAX_CPU_NUMBER=12',
        '-DMAX_PARALLEL_NUMBER=1',
        '-DBUILD_SINGLE=1',
        '-DBUILD_DOUBLE=1',
        '-DBUILD_COMPLEX=1',
        '-DBUILD_COMPLEX16=1',
        '-UASMNAME',
        '-UASMFNAME',
        '-UNAME',
        '-UCNAME',
        '-UCHAR_NAME',
        '-UCHAR_CNAME',
        '-DASMNAME=cblas_drot',
        '-DASMFNAME=cblas_drot_',
        '-DNAME=cblas_drot_',
        '-DCNAME=cblas_drot',
        '-DCHAR_NAME="cblas_drot_"',
        '-DCHAR_CNAME="cblas_drot"',
        '-DNO_AFFINITY',
        '-DDOUBLE',
        '-UCOMPLEX',
        '-DCBLAS',],
)

_interface = static_library('_interface',
                            c_args: _cargs,
                            link_whole: _static_libs)
