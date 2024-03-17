# Ordered As per https://netlib.org/blas/blasqr.pdf
# NOTE: xROTG xROTMG xROTM have no kernels?
_kinds = []
_rk = ['s', 'd']
_ck = ['c', 'z']
_std_dk = _rk + _ck # Standard precisions
# TODO: Actually test and set this
if true
  fma3_flag = '-mfma'
endif
# TODO: This is currently following x86_64 generic for src and dir, but it needs
# to diversify
# NOTE: The def and undefs are from Makefile.L1
# Construct all PRECsymbKERNEL from src and dir via files(dir + src)
# For the precs array, the following mapping is used for c_args:
# addl --> passed AS IS
base_kops = [
  # Level 1 BLAS
  {
    'base': 'rot',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'rot.c',
        'addl': [fma3_flag],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'rot.c',
        'addl': [fma3_flag],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
      },
      'cs' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
      },
      'zd' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
      },
      'xq' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
      },
    },
  },
  {'base': 'swap',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'swap.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'swap.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
      },
    },
  },
  {'base': 'scal',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'scal.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'scal.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
      },
    },
  },
  {'base': 'copy',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'copy.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'copy.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
      },
    },
  },
  {'base': 'axpy',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'axpy.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'axpy.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
      },
    },
  },
  {'base': 'dot',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'dot.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'dot.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
      },
    },
  },
  # xDOTU xDOTC xxDOT aren't present
  {'base': 'nrm2',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'nrm2.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'nrm2.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
      },
    },
  },
  {'base': 'asum',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'asum.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'asum.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
      },
    },
  },
  {'base': 'amax',
    'precs': {
      's' : {
        'dir': 'arm',
        'kernel': 'amax.c',
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'amax.c',
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
      },
    },
  },
]

kernel_confs = []

# For the precs array, the following mapping is used for c_args:
# undef --> -Uwhatever
# def --> -Dwhatever
foreach _kop : base_kops
  foreach pkey, pval : _kop['precs']
    kcfg = {
      'name': pkey + _kop['base'] + '_k',
      # TODO: This should probably be in files() after this is ready
      # dictionaries with files can't be printed with message
      'src': pval['dir'] + '/' + pval['kernel'],
    }
    if pval.has_key('addl')
      kcfg += {'addl': pval['addl']}
    endif
    if 's' == pkey
      kcfg += {'undef': ['COMPLEX', 'DOUBLE']}
    endif
    if 'd' == pkey
      kcfg += {
        'undef': ['COMPLEX'],
        'def': ['DOUBLE'],
      }
    endif
    if 'q' == pkey
      kcfg += {
        'undef': ['COMPLEX'],
        'def': ['DXDOUBLE'],
      }
    endif
    # TODO: rot is nonstandard in this instance, taking cs zd xq
    # The others (scal, swap) take c z x, so we need to use in
    if 'c' in pkey
      kcfg += {
        'undef': ['DOUBLE'],
        'def': ['COMPLEX'],
      }
    endif
    if 'z' in pkey
      kcfg += {
        'def': ['COMPLEX', 'DOUBLE'],
      }
    endif
    if 'x' in pkey
      kcfg += {
        'def': ['COMPLEX', 'DXDOUBLE'],
      }
    endif
    message(kcfg)
  endforeach
endforeach

# foreach root : blas1_roots
#   fname = root + '.S'
#   defs = []
#   foreach prec : real_kinds
#     name = prec + fname + '_k'
#     kernel_confs += {'defs': defs, 'name': name, 'src': fname}
#   endforeach
# endforeach

# _static_libs = []
# foreach conf: kernel_confs
#     _static_libs += static_library(
#         conf['name'],
#         conf['src'],
#         include_directories:  _inc,
#         c_args: conf['defs'],
#     )
# endforeach
