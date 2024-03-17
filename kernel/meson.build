# Ordered As per https://netlib.org/blas/blasqr.pdf
# NOTE: xROTG xROTMG xROTM have no kernels?
_kinds = []
_rk = ['s', 'd']
_ck = ['c', 'z']
_std_dk = _rk + _ck # Standard precisions
# TODO: Actually test and set this
if fma3
  fma3_flag = '-mfma'
endif
# TODO: This is currently following x86_64 generic for src and dir, but it needs
# to diversify
# NOTE: The def and undefs are from Makefile.L1
# Construct all PRECsymbKERNEL from src and dir via files(dir + src)
# For the precs array, the following mapping is used for c_args:
# undef --> -Uwhatever
# def --> -Dwhatever
# addl --> passed AS IS
kops_array = [
  # Level 1 BLAS
  {
    'base': 'rot',
    # TODO: A lot of these flags are basically repeating..
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'rot.c',
        'undef': ['COMPLEX', 'DOUBLE'],
        'addl': [fma3_flag],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'rot.c',
        'undef': ['COMPLEX'],
        'def': ['DOUBLE'],
        'addl': [fma3_flag],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'undef': ['COMPLEX'],
        'def': ['DXDOUBLE'],
      },
      'cs' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'undef': ['DOUBLE'],
        'def': ['COMPLEX'],
      },
      'zd' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'def': ['COMPLEX', 'DOUBLE'],
      },
      'xq' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'def': ['COMPLEX', 'DXDOUBLE'],
      },
    },
  },
  {'symb': 'swap',
   'src': {'real': 'swap', 'cmplx': 'zswap', 'ext': '.c'},
   'dir': 'arm', 'precs': _std_dk },
  {'symb': 'scal',
   'src': {'real': 'scal', 'cmplx': 'zscal', 'ext': '.c'},
   'dir': 'arm', 'precs': _std_dk },
  {'symb': 'copy',
   'src': {'real': 'copy', 'cmplx': 'zcopy', 'ext': '.c'},
   'dir': 'arm', 'precs': _std_dk },
  {'symb': 'axpy',
   'src': {'real': 'axpy', 'cmplx': 'zaxpy', 'ext': '.c'},
   'dir': 'arm', 'precs': _std_dk },
  {'symb': 'dot',
   'src': {'real': 'dot', 'cmplx': 'zdot', 'ext': '.c'},
   'dir': 'arm', 'precs': _std_dk },
]

kernel_confs = []

foreach root : blas1_roots
  fname = root + '.S'
  defs = []
  foreach prec : real_kinds
    name = prec + fname + '_k'
    kernel_confs += {'defs': defs, 'name': name, 'src': fname}
  endforeach
endforeach

_static_libs = []
foreach conf: kernel_confs
    _static_libs += static_library(
        conf['name'],
        conf['src'],
        include_directories:  _inc,
        c_args: conf['defs'],
    )
endforeach
