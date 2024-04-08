# Ordered As per https://netlib.org/blas/blasqr.pdf
# NOTE: xROTG xROTMG xROTM have no kernels?
# TODO: Actually test and set this
if true
  fma3_flag = '-mfma'
endif
# TODO: This is currently following x86_64 generic for src and dir, but it needs
# to diversify
# NOTE: The def and undefs are from Makefile.L1
# Construct all ModesymbKERNEL from src and dir via files(dir + src)
# For the modes array, the following mapping is used for c_args:
# addl --> passed AS IS
base_kops = [
  # Level 1 BLAS
  {
    'base': 'rot',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'rot.c',
        'addl': [fma3_flag],
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'rot.c',
        'addl': [fma3_flag],
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'exts': ['_k'],
      },
      'cs' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'exts': ['_k'],
      },
      'zd' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'exts': ['_k'],
      },
      'xq' : {
        'dir': 'arm',
        'kernel': 'zrot.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'swap',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'swap.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'swap.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zswap.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'scal',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'scal.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'scal.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zscal.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'copy',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'copy.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'copy.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zcopy.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'axpy',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'axpy.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'axpy.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zaxpy.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'dot',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'dot.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'dot.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zdot.c',
        'exts': ['_k'],
      },
    },
  },
  # xDOTU xDOTC xxDOT aren't present
  {'base': 'nrm2',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'nrm2.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'nrm2.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'znrm2.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'asum',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'asum.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'asum.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zasum.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'amax',
    'modes': {
      's' : {
        'dir': 'arm',
        'kernel': 'amax.c',
        'exts': ['_k'],
      },
      'd' : {
        'dir': 'arm',
        'kernel': 'amax.c',
        'exts': ['_k'],
      },
      'q' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
        'exts': ['_k'],
      },
      'c' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
        'exts': ['_k'],
      },
      'z' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
        'exts': ['_k'],
      },
      'x' : {
        'dir': 'arm',
        'kernel': 'zamax.c',
        'exts': ['_k'],
      },
    },
  },
  {'base': 'axpby',
   'modes': {'s': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
             'd': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
             'c': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
             'z': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
            },
  },
  # Level 2 BLAS
  # There are additional sources so now we have srcs
  # Ordered as per KERNEL.generic and Makefile.L2
  # exts are used to find the flags for each modality
  # ext is attached to base (only useful from Level 2)
  # i.e. baseext (e.g., gemv_n, gemv_t, cgeru_k, cgerc_k)
  {'base': 'gemv',
   'modes': {'s': {'dir': 'arm', 'kernel': 'gemv_n.c', 'exts': ['_n', '_t']},
             'd': {'dir': 'arm', 'kernel': 'gemv_n.c',
                   # TODO: _t should take fmaflag, but then we need a dictionary..
                   # From Makefile.L2
                   'exts': ['_n', '_t']},
             'q': {'dir': 'arm', 'kernel': 'gemv_n.c', 'exts': ['_n', '_t']},
             'c': {'dir': 'arm', 'kernel': 'zgemv_n.c',
                   'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
             'z': {'dir': 'arm', 'kernel': 'zgemv_n.c',
                   'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
             'x': {'dir': 'arm', 'kernel': 'zgemv_n.c',
                   'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
            },
  },
  {'base': 'symv',
   'modes': {'s': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
             'd': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
             'q': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
             'c': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
             'z': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
             'x': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
            },
  },
  {'base': 'ger',
   'modes': {'s': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']},
             'd': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']},
             'q': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']}
            },
  },
  {'base': 'geru',
   'modes': {'c': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
             'z': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
             'x': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']}
            }
  },
  {'base': 'gerc',
   'modes': {'c': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
             'z': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
             'x': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
            },
  },
  {'base': 'hemv',
   'modes': {'c': {'dir': 'generic',
                   'kernel': 'zhemv_k.c',
                   'exts': ['_U', '_L', '_V', '_M']},
             'z': {'dir': 'generic',
                   'kernel': 'zhemv_k.c',
                   'exts': ['_U', '_L', '_V', '_M']},
             'x': {'dir': 'generic',
                   'kernel': 'zhemv_k.c',
                   'exts': ['_U', '_L', '_V', '_M']},
            },
  },
  {'base': 'bgemv',
   'modes': {'s': {'dir': 'x86_64',
                   'kernel': 'sbgemv_n.c',
                   'exts': ['_n', '_t']}
            },
  },
]

# kernel_confs = []

# # For the modes array, the following mapping is used for c_args:
# # undef --> -Uwhatever
# # def --> -Dwhatever
# foreach _kop : base_kops
#   foreach pkey, pval : _kop['modes']
#     # TODO: This srcs should be in files() after this is ready
#     # dictionaries with files can't be printed with message
#     kcfg = {
#       'src': pval['dir'] + '/' + pval['kernel'],
#     }
#     if pval.has_key('addl')
#       kcfg += {'addl': pval['addl']}
#     endif
#     if 's' == pkey
#       foreach ext : pval['exts']
#         kcfg += {'name': pkey + _kop['base'] + ext}
#         # Level 1
#         if '_k' == ext
#           kcfg += {
#             'undef': ['COMPLEX', 'DOUBLE']
#           }
#         endif
#         # Level 2
#         if '_n' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'COMPLEX', 'TRANS'],
#           }
#         endif
#         if '_t' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'COMPLEX'],
#             'def': ['TRANS']
#           }
#         endif
#         if '_U' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'COMPLEX', 'LOWER'],
#           }
#         endif
#         if '_L' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'COMPLEX'],
#             'def': ['LOWER']
#           }
#         endif
#       endforeach
#     endif
#     if 'd' == pkey
#       foreach ext : pval['exts']
#         kcfg += {'name': pkey + _kop['base'] + ext}
#         # Level 1
#         if '_k' == ext
#           kcfg += {
#             'undef': ['COMPLEX'],
#             'def': ['DOUBLE'],
#           }
#         endif
#         # Level 2
#         if '_n' == ext
#           kcfg += {
#             'undef': ['COMPLEX', 'TRANS'],
#             'def': ['DOUBLE']
#           }
#         endif
#         if '_t' == ext
#           kcfg += {
#             'undef': ['COMPLEX'],
#             'def': ['DOUBLE', 'TRANS']
#           }
#         endif
#         if '_U' == ext
#           kcfg += {
#             'undef': ['COMPLEX', 'LOWER'],
#             'def': ['DOUBLE']
#           }
#         endif
#         if '_L' == ext
#           kcfg += {
#             'undef': ['COMPLEX'],
#             'def': ['DOUBLE', 'LOWER']
#           }
#         endif
#       endforeach
#     endif
#     if 'q' == pkey
#       foreach ext : pval['exts']
#         kcfg += {'name': pkey + _kop['base'] + ext}
#         # Level 1
#         if '_k' == ext
#           kcfg += {
#             'undef': ['COMPLEX'],
#             'def': ['XDOUBLE'],
#           }
#         endif
#         # Level 2
#         if '_n' == ext
#           kcfg += {
#             'undef': ['COMPLEX', 'TRANS'],
#             'def': ['XDOUBLE']
#           }
#         endif
#         if '_t' == ext
#           kcfg += {
#             'undef': ['COMPLEX'],
#             'def': ['XDOUBLE', 'TRANS']
#           }
#         endif
#         if '_U' == ext
#           kcfg += {
#             'undef': ['COMPLEX', 'LOWER'],
#             'def': ['XDOUBLE']
#           }
#         endif
#         if '_L' == ext
#           kcfg += {
#             'undef': ['COMPLEX'],
#             'def': ['XDOUBLE', 'LOWER']
#           }
#         endif
#       endforeach
#     endif
#     # TODO: rot is nonstandard in this instance, taking cs zd xq
#     # The others (scal, swap) take c z x, so we need to use in
#     if 'c' in pkey
#       foreach ext : pval['exts']
#         kcfg += {'name': pkey + _kop['base'] + ext}
#         # Level 1
#         if '_k' == ext
#           kcfg += {
#             'undef': ['DOUBLE'],
#             'def': ['COMPLEX'],
#           }
#         endif
#         # Level 2
#         if '_n' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'TRANS', 'CONJ', 'XCONJ'],
#             'def': ['COMPLEX']
#           }
#         endif
#         if '_t' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'CONJ', 'XCONJ'],
#             'def': ['COMPLEX', 'TRANS']
#           }
#         endif
#         if '_r' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'TRANS', 'XCONJ'],
#             'def': ['COMPLEX', 'CONJ']
#           }
#         endif
#         if '_c' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'XCONJ'],
#             'def': ['COMPLEX', 'TRANS', 'CONJ']
#           }
#         endif
#         if '_o' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'TRANS', 'CONJ'],
#             'def': ['COMPLEX', 'XCONJ']
#           }
#         endif
#         if '_u' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'CONJ'],
#             'def': ['COMPLEX', 'TRANS', 'XCONJ']
#           }
#         endif
#         if '_s' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'TRANS'],
#             'def': ['COMPLEX', 'CONJ', 'XCONJ']
#           }
#         endif
#         if '_d' == ext
#           kcfg += {
#             'undef': ['DOUBLE'],
#             'def': ['COMPLEX', 'TRANS', 'CONJ', 'XCONJ']
#           }
#         endif
#         if '_U' == ext
#           kcfg += {
#             'undef': ['DOUBLE', 'LOWER'],
#             'def': ['COMPLEX']
#           }
#         endif
#         if '_L' == ext
#           kcfg += {
#             'undef': ['DOUBLE'],
#             'def': ['COMPLEX', 'LOWER']
#           }
#         endif
#       endforeach
#     endif
#     if 'z' in pkey
#       foreach ext : pval['exts']
#         kcfg += {'name': pkey + _kop['base'] + ext}
#         # Level 1
#         if '_k' == ext
#           kcfg += {
#             'def': ['COMPLEX', 'DOUBLE'],
#           }
#         endif
#         # Level 2
#         if '_n' == ext
#           kcfg += {
#             'undef': ['CONJ', 'XCONJ', 'TRANS'],
#             'def': ['COMPLEX', 'DOUBLE'],
#           }
#         endif
#         if '_t' == ext
#           kcfg += {
#             'undef': ['CONJ', 'XCONJ'],
#             'def': ['COMPLEX', 'DOUBLE', 'TRANS']
#           }
#         endif
#         if '_r' == ext
#           kcfg += {
#             'undef': ['TRANS', 'XCONJ'],
#             'def': ['COMPLEX', 'DOUBLE', 'CONJ']
#           }
#         endif
#         if '_c' == ext
#           kcfg += {
#             'undef': ['XCONJ'],
#             'def': ['COMPLEX', 'DOUBLE', 'TRANS', 'CONJ']
#           }
#         endif
#         if '_o' == ext
#           kcfg += {
#             'undef': ['CONJ', 'TRANS'],
#             'def': ['COMPLEX', 'DOUBLE', 'XCONJ']
#           }
#         endif
#         if '_u' == ext
#           kcfg += {
#             'undef': ['CONJ'],
#             'def': ['COMPLEX', 'DOUBLE', 'TRANS', 'XCONJ']
#           }
#         endif
#         if '_s' == ext
#           kcfg += {
#             'undef': ['TRANS'],
#             'def': ['COMPLEX', 'DOUBLE', 'CONJ', 'XCONJ']
#           }
#         endif
#         if '_d' == ext
#           kcfg += {
#             'def': ['COMPLEX', 'DOUBLE', 'TRANS', 'CONJ', 'XCONJ']
#           }
#         endif
#         if '_U' == ext
#           kcfg += {
#             'undef': ['LOWER'],
#             'def': ['COMPLEX', 'DOUBLE']
#           }
#         endif
#         if '_L' == ext
#           kcfg += {
#             'def': ['COMPLEX', 'DOUBLE', 'LOWER']
#           }
#         endif
#       endforeach
#     endif
#     if 'x' in pkey
#       foreach ext : pval['exts']
#         kcfg += {'name': pkey + _kop['base'] + ext}
#         # Level 1
#         if '_k' == ext
#           kcfg += {
#             'def': ['COMPLEX', 'XDOUBLE'],
#           }
#         endif
#         # Level 2
#         if '_n' == ext
#           kcfg += {
#             'undef': ['XCONJ', 'CONJ', 'TRANS'],
#             'def': ['XDOUBLE', 'COMPLEX'],
#           }
#         endif
#         if '_t' == ext
#           kcfg += {
#             'undef': ['CONJ', 'XCONJ'],
#             'def': ['XDOUBLE', 'COMPLEX', 'TRANS']
#           }
#         endif
#         if '_r' == ext
#           kcfg += {
#             'undef': ['TRANS', 'XCONJ'],
#             'def': ['XDOUBLE', 'COMPLEX', 'CONJ']
#           }
#         endif
#         if '_c' == ext
#           kcfg += {
#             'undef': ['XCONJ'],
#             'def': ['XDOUBLE', 'COMPLEX', 'TRANS', 'CONJ']
#           }
#         endif
#         if '_o' == ext
#           kcfg += {
#             'undef': ['TRANS', 'CONJ'],
#             'def': ['XDOUBLE', 'COMPLEX', 'XCONJ']
#           }
#         endif
#         if '_u' == ext
#           kcfg += {
#             'undef': ['CONJ'],
#             'def': ['XDOUBLE', 'COMPLEX', 'TRANS', 'XCONJ']
#           }
#         endif
#         if '_s' == ext
#           kcfg += {
#             'undef': ['TRANS'],
#             'def': ['XDOUBLE', 'COMPLEX', 'CONJ', 'XCONJ']
#           }
#         endif
#         if '_d' == ext
#           kcfg += {
#             'def': ['XDOUBLE', 'COMPLEX', 'TRANS', 'CONJ', 'XCONJ']
#           }
#         endif
#         if '_U' == ext
#           kcfg += {
#             'undef': ['LOWER'],
#             'def': ['COMPLEX', 'XDOUBLE']
#           }
#         endif
#         if '_L' == ext
#           kcfg += {
#             'def': ['COMPLEX', 'XDOUBLE', 'LOWER']
#           }
#         endif
#       endforeach
#     endif
#     message(kcfg)
#   endforeach
# endforeach

# # foreach root : blas1_roots
# #   fname = root + '.S'
# #   defs = []
# #   foreach prec : real_kinds
# #     name = prec + fname + '_k'
# #     kernel_confs += {'defs': defs, 'name': name, 'src': fname}
# #   endforeach
# # endforeach

# # _static_libs = []
# # foreach conf: kernel_confs
# #     _static_libs += static_library(
# #         conf['name'],
# #         conf['src'],
# #         include_directories:  _inc,
# #         c_args: conf['defs'],
# #     )
# # endforeach
# base_kops = [
#   {'base': 'axpby',
#    'modes': {'s': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
#              'd': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
#              'c': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
#              'z': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
#             },
#    }
# ]

_configs = [
  # {'name': 'saxpby_k', 'undef': ['COMPLEX', 'DOUBLE'], 'def': [], 'addl': []},
  # {'name': 'daxpby_k', 'undef': ['COMPLEX'], 'def': ['DOUBLE'], 'addl': []},
  # {'name': 'caxpby_k', 'undef': ['CONJ', 'DOUBLE'], 'def': ['COMPLEX'], 'addl': []},
  # {'name': 'zaxpby_k', 'undef': ['CONJ'], 'def': ['COMPLEX', 'DOUBLE'], 'addl': []},
  # {'name': 'sswap_k', 'undef': ['COMPLEX', 'DOUBLE'], 'def': [], 'addl': []},
  # {'name': 'dswap_k', 'undef': ['COMPLEX'], 'def': ['DOUBLE'], 'addl': []},
  # {'name': 'qswap_k', 'undef': ['COMPLEX'], 'def': ['XDOUBLE'], 'addl': []},
  # {'name': 'cswap_k', 'undef': ['DOUBLE'], 'def': ['COMPLEX'], 'addl': []},
  # {'name': 'zswap_k', 'undef': [], 'def': ['COMPLEX', 'DOUBLE'], 'addl': []},
  # {'name': 'xswap_k', 'undef': [], 'def': ['COMPLEX', 'XDOUBLE'], 'addl': []},
  # {'name': 'srot_k',
  #  'undef': ['COMPLEX', 'COMPLEX', 'DOUBLE'],
  #  'def': [],
  #  'addl': [fma3_flag]},
  {'name': 'drot_k',
   'undef': ['COMPLEX', 'COMPLEX'],
   'def': ['DOUBLE'],
   'addl': [fma3_flag,
        '-m64',
        '-DSMP_SERVER',
        '-DNO_WARMUP',
        '-DMAX_CPU_NUMBER=12',
        '-DMAX_PARALLEL_NUMBER=1',
        '-DBUILD_SINGLE=1',
        '-DBUILD_DOUBLE=1',
        '-DBUILD_COMPLEX=1',
        '-DBUILD_COMPLEX16=1',
        '-msse3',
        '-mssse3',
        '-msse4.1',
        '-mavx',
        '-mavx2',
        '-UASMNAME',
        '-UASMFNAME',
        '-UNAME',
        '-UCNAME',
        '-UCHAR_NAME',
        '-UCHAR_CNAME',
        '-DASMNAME=drot_k',
        '-DASMFNAME=drot_k_',
        '-DNAME=drot_k_',
        '-DCNAME=drot_k',
        '-DCHAR_NAME="drot_k_"',
        '-DCHAR_CNAME="drot_k"',
        '-DNO_AFFINITY',
        '-DDOUBLE',
        '-UCOMPLEX',
        '-DCBLAS',
           ]},
  # {'name': 'qrot_k',
  #  'undef': ['COMPLEX', 'COMPLEX'],
  #  'def': ['XDOUBLE'],
  #  'addl': []},
  # {'name': 'csrot_k',
  #  'undef': ['DOUBLE'],
  #  'def': ['COMPLEX', 'COMPLEX'],
  #  'addl': []},
  # {'name': 'zdrot_k',
  #  'undef': [],
  #  'def': ['COMPLEX', 'COMPLEX', 'DOUBLE'],
  #  'addl': []},
  # {'name': 'xqrot_k',
  #  'undef': [],
  #  'def': ['COMPLEX', 'COMPLEX', 'XDOUBLE'],
  #  'addl': []}
]

kernel_confs = []

foreach _kop : base_kops
  base = _kop['base']
  modes = _kop['modes']
  foreach mode, details : modes
    dir = details['dir']
    kernel = details['kernel']
    exts = details['exts']
    foreach ext : exts
      src = join_paths(dir, kernel)
      foreach cfg : _configs
        if mode + base + ext == cfg['name']
          defs = cfg.get('def', [])
          undefs = cfg.get('undef', [])
          addl = cfg.get('addl', [])
          _cargs = []
          foreach _d : defs
            _cargs += ('-D' + _d)
          endforeach
          foreach _u : undefs
            _cargs += ('-U' + _u)
          endforeach
          foreach _addl : addl
            _cargs += addl
          endforeach
          current_def = {'c_args': _cargs, 'name': cfg['name'], 'src': src}
          message(current_def)
          kernel_confs += current_def
        endif
      endforeach
    endforeach
  endforeach
endforeach

_kern_libs = []
foreach conf: kernel_confs
    _kern_libs += static_library(
        conf['name'],
        conf['src'],
        include_directories:  _inc,
        c_args: conf['c_args'],
    )
endforeach

_kern = static_library('_kern',
                link_whole: _kern_libs)
