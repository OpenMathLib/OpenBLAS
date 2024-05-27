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
  { 'base': '?rot',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'rot.c', 'addl': [fma3_flag]}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'rot.c', 'addl': [fma3_flag]}}},
        'q': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zrot.c'}}},
        'cs': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zrot.c'}}},
        'zd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zrot.c'}}},
        'xq': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zrot.S'}}},
    },
  },
  { 'base': '?swap',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'swap.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'swap.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zswap.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zswap.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'swap.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zswap.S'}}},
    },
  },
  { 'base': '?scal',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'scal.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'scal.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zscal.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zscal.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'scal.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zscal.S'}}},
    },
  },
  { 'base': '?copy',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'copy.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'copy.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zcopy.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zcopy.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'copy.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zcopy.S'}}},
    },
  },
  { 'base': '?axpy',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'axpy.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'axpy.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zaxpy.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zaxpy.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'axpy.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zaxpy.S'}}},
    },
  },
  { 'base': '?dot',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'dot.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'dot.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zdot.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zdot.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'dot.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zdot.S'}}},
    },
  },
  # TODO(rg): Add dsdotkernel conditionals
  # xDOTU xDOTC xxDOT aren't present
  { 'base': '?nrm2',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'nrm2.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'nrm2.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'znrm2.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'znrm2.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'nrm2.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'znrm2.S'}}},
    },
  },
  { 'base': '?asum',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'asum.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'asum.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zasum.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zasum.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'asum.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zasum.S'}}},
    },
  },
  { 'base': '?amax',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'amax.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'amax.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamax.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamax.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zamax.S'}}},
    },
  },
  { 'base': '?sum',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'sum.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'sum.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zsum.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zsum.c'}}},
    },
  },
  { 'base': '?amin',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'amin.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'amin.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamin.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamin.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'izamax.S'}}},
    },
  },
  { 'base': 'i?amax',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'iamax.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'iamax.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'izamax.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'izamax.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
        'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'izamax.S'}}},
    },
  },
  { 'base': 'i?amin',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'iamin.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'iamin.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'izamin.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'izamin.c'}}},
    },
  },
  { 'base': 'i?max',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'imax.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'imax.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
    },
  },
  { 'base': 'i?min',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'imin.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'imin.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
    },
  },
  { 'base': '?max',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'max.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'max.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax.S'}}},
    },
  },
  { 'base': '?min',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'min.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'min.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax.S'}}},
    },
  },
  { 'base': '?axpby',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'axpby.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'axpby.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zaxpby.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zaxpby.c'}}},
    },
  },
  # Level 2 BLAS
  # There are additional sources so now we have srcs
  # Ordered as per KERNEL.generic and Makefile.L2
  # exts are used to find the flags for each modality
  # ext is attached to base (only useful from Level 2)
  # i.e. baseext (e.g., gemv_n, gemv_t, cgeru_k, cgerc_k)
  { 'base': '?gemv',
    'modes': {
        's': {
            'exts': {
              # TODO(rg): Where are these coming from??
              # Most of these have both generic defines and also per-folder defines..
              # Makefile lists sgemv_n_4.c as the source, though there is a sgemv_n.c
                '_n': {'dir': 'arm', 'kernel': 'gemv_n.c'},
                '_t': {'dir': 'arm', 'kernel': 'gemv_t.c'},
            }
        },
        'd': {
            'exts': {
                '_n': {'dir': 'arm', 'kernel': 'gemv_n.c'},
                '_t': {'dir': 'arm', 'kernel': 'gemv_t.c'},
            }
        },
        'q': {
            'exts': {
                '_n': {'dir': 'x86_64', 'kernel': 'qgemv_n.S'},
                '_t': {'dir': 'x86_64', 'kernel': 'qgemv_t.S'},
            }
        },
        'c': {
            'exts': {
                '_n': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_t': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_r': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_c': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_o': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_u': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_s': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_d': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
            }
        },
        'z': {
            'exts': {
                '_n': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_t': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_r': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_c': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_o': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_u': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_s': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_d': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
            }
        },
        'x': {
            'exts': {
                '_n': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_t': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_r': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_c': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_o': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_u': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_s': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
                '_d': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
            }
        },
    },
  },
  { 'base': '?symv',
    'modes': {
        's': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'symv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'symv_k.c'},
            }
        },
        'd': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'symv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'symv_k.c'},
            }
        },
        'c': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
            }
        },
        'z': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
            }
        },
        'q': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'symv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'symv_k.c'},
            }
        },
        'x': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
            }
        },
    },
  },
  { 'base': '?lsame',
    'modes': {
        '': {'exts': {'': {'dir': 'generic', 'kernel': 'lsame.c'}}},
    },
  },
  { 'base': '?cabs',
    'modes': {
        's': {'exts': {'': {'dir': 'generic', 'kernel': 'cabs.c'}}},
        'd': {'exts': {'': {'dir': 'generic', 'kernel': 'cabs.c'}}},
        'q': {'exts': {'': {'dir': 'generic', 'kernel': 'cabs.c'}}},
    },
  },
  # TODO(rg): There's CGEMM3MKERNEL, ZGEMM3MKERNEL in kernel.generic
  { 'base': '?ger',
    'modes': {
        's': {'exts': {'_k': {'dir': 'generic', 'kernel': 'ger.c'}}},
        'd': {'exts': {'_k': {'dir': 'generic', 'kernel': 'ger.c'}}},
        'q': {'exts': {'_k': {'dir': 'generic', 'kernel': 'ger.c'}}},
    },
  },
  { 'base': '?geru',
    'modes': {
        'c': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        'z': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        'x': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
    },
  },
  { 'base': '?gerc',
    'modes': {
        'c': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        'z': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        'x': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
    },
  },
  { 'base': '?hemv',
    'modes': {
        'c': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_V': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_M': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
            }
        },
        'z': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_V': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_M': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
            }
        },
        'x': {
            'exts': {
                '_U': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_L': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_V': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
                '_M': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
            }
        },
    },
  },
  { 'base': '?bgemv',
    'modes': {
        's': {
            'exts': {
                '_n': {'dir': 'x86_64', 'kernel': 'sbgemv_n.c'},
                '_t': {'dir': 'x86_64', 'kernel': 'sbgemv_n.c'},
            }
        }
    },
  },
  # Level 3 symbols
  { 'base': '?gemm_kernel',
    'modes': {
        's': {'exts': {'': {'dir': 'generic', 'kernel': 'gemmkernel_2x2.c'}}},
        'd': {'exts': {'': {'dir': 'generic', 'kernel': 'gemmkernel_2x2.c'}}},
        'c': {
            'exts': {
                '_n': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DNN']},
                '_l': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DCN']},
                # TODO(rg): What about _r conditionals? Makefile.L3:2969
                '_r': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DNC']},
                '_b': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DCC']},
            }
        },
        'z': {
            'exts': {
                '_n': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DNN']},
                '_l': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DCN']},
                '_r': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DNC']},
                '_b': {'dir': 'generic', 'kernel': 'zgemmkernel_2x2.c', 'addl': ['-DCC']},
            }
        }
        # 'q': {'exts': {'': {'dir': 'generic', 'kernel': 'gemm_beta.c'}}},
        # 'x': {'exts': {'': {'dir': 'generic', 'kernel': 'zgemm_beta.c'}}},
    },
  },
  { 'base': '?trmm_kernel',
    'modes': {
        's': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c'},
                '_LT': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c', 'addl': ['-DLEFT', '-DTRANSA']},
                '_RN': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c'},
                '_RT': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c'},
            }
        },
        'd': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c'},
                '_LT': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c', 'addl': ['-DLEFT', '-DTRANSA']},
                '_RN': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c'},
                '_RT': {'dir': 'generic', 'kernel': 'trmmkernel_2x2.c'},
            }
        },
        'c': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-UCONJ', '-DNN']},
                '_LT': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DLEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_LR': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DCN']},
                '_LC': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DCN']},
                '_RN': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-UCONJ', '-DNN']},
                '_RT': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-ULEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_RR': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DNC']},
                '_RC': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DCN']},
            }
        },
        'z': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-UCONJ', '-DNN']},
                '_LT': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DLEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_LR': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DCN']},
                '_LC': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DCN']},
                '_RN': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-UCONJ', '-DNN']},
                '_RT': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-ULEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_RR': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DNC']},
                '_RC': {'dir': 'generic', 'kernel': 'ztrmmkernel_2x2.c',
                        'addl': ['-DCONJ', '-DCN']},
            },
        },
    },
  },
  { 'base': '?trsm_kernel',
    'modes': {
        's': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'trsm_kernel_LN.c',
                        'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'dir': 'generic', 'kernel': 'trsm_kernel_LT.c',
                        'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_RN': {'dir': 'generic', 'kernel': 'trsm_kernel_RN.c',
                        'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'dir': 'generic', 'kernel': 'trsm_kernel_RT.c',
                        'addl': ['-DRT', '-UUPPER', '-UCONJ']},
            },
        },
        'd': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'trsm_kernel_LN.c',
                        'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'dir': 'generic', 'kernel': 'trsm_kernel_LT.c',
                        'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_RN': {'dir': 'generic', 'kernel': 'trsm_kernel_RN.c',
                        'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'dir': 'generic', 'kernel': 'trsm_kernel_RT.c',
                        'addl': ['-DRT', '-UUPPER', '-UCONJ']},
            },
        },
        'c': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'trsm_kernel_LN.c',
                        'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'dir': 'generic', 'kernel': 'trsm_kernel_LT.c',
                        'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_LR': {'dir': 'generic', 'kernel': 'trsm_kernel_LN.c',
                        'addl': ['-DLN', '-DUPPER', '-DCONJ']},
                '_LC': {'dir': 'generic', 'kernel': 'trsm_kernel_LT.c',
                        'addl': ['-DLT', '-UUPPER', '-DCONJ']},
                '_RN': {'dir': 'generic', 'kernel': 'trsm_kernel_RN.c',
                        'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'dir': 'generic', 'kernel': 'trsm_kernel_RT.c',
                        'addl': ['-DRT', '-UUPPER', '-UCONJ']},
                '_RR': {'dir': 'generic', 'kernel': 'trsm_kernel_RN.c',
                        'addl': ['-DRN', '-DUPPER', '-DCONJ']},
                '_RC': {'dir': 'generic', 'kernel': 'trsm_kernel_RT.c',
                        'addl': ['-DRT', '-UUPPER', '-DCONJ']},
            },
        },
    },
  },
  { 'base': '?gemm3m_kernel',
    'modes': {
        'c': {'exts': {'': {'dir': 'x86_64', 'kernel': 'zgemm3m_kernel_8x4_sse3.S',
                           'addl': ['-DNN']}}},
        'z': {'exts': {'': {'dir': 'x86_64', 'kernel': 'zgemm3m_kernel_8x4_sse3.S',
                           'addl': ['-DNN']}}},
    },
  },
  { 'base': '?gemm',
    'modes': {
        's': {'exts': {
            '_beta': {'dir': 'generic', 'kernel': 'gemm_beta.c'},
            # TODO(rg): the _NUM prefixes are arch dependent
            '_incopy': {'dir': 'generic', 'kernel': 'gemm_ncopy_8.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_8.c'},
            '_oncopy': {'dir': 'generic', 'kernel': 'gemm_ncopy_2.c'},
            '_otcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_4.c'},
            # TODO(rg): direct and direct_performant are built only conditionally
            '_direct': {'dir': 'x86_64', 'kernel': 'sgemm_direct_skylakex.c'},
            '_direct_performant': {'dir': 'x86_64', 'kernel': 'sgemm_direct_performant.c'},
        }},
        'd': {'exts': {
            '_beta': {'dir': 'generic', 'kernel': 'gemm_beta.c'},
            '_incopy': {'dir': 'generic', 'kernel': 'gemm_ncopy_8.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_8.c'},
            '_oncopy': {'dir': 'generic', 'kernel': 'gemm_ncopy_2.c'},
            '_otcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_4.c'},
        }},
        'c': {'exts': {
            '_beta': {'dir': 'generic', 'kernel': 'zgemm_beta.c'},
            '_incopy': {'dir': 'generic', 'kernel': 'zgemm_ncopy_8.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'zgemm_tcopy_8.c'},
            '_oncopy': {'dir': 'generic', 'kernel': 'zgemm_ncopy_2.c'},
            '_otcopy': {'dir': 'generic', 'kernel': 'zgemm_tcopy_2.c'},
        }},
        'z': {'exts': {
            '_beta': {'dir': 'generic', 'kernel': 'zgemm_beta.c'},
            '_incopy': {'dir': 'generic', 'kernel': 'zgemm_ncopy_8.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'zgemm_tcopy_8.c'},
            '_oncopy': {'dir': 'generic', 'kernel': 'zgemm_ncopy_2.c'},
            '_otcopy': {'dir': 'generic', 'kernel': 'zgemm_tcopy_2.c'},
        }},
    },
  },
  { 'base': '?trmm',
    'modes': {
        's': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'dir': 'generic', 'kernel': 'trmm_uncopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'dir': 'generic', 'kernel': 'trmm_uncopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'dir': 'generic', 'kernel': 'trmm_lncopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'dir': 'generic', 'kernel': 'trmm_lncopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'dir': 'generic', 'kernel': 'trmm_utcopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'dir': 'generic', 'kernel': 'trmm_utcopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'dir': 'generic', 'kernel': 'trmm_ltcopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'dir': 'generic', 'kernel': 'trmm_ltcopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'dir': 'generic', 'kernel': 'trmm_uncopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'dir': 'generic', 'kernel': 'trmm_uncopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'dir': 'generic', 'kernel': 'trmm_lncopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'dir': 'generic', 'kernel': 'trmm_lncopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'dir': 'generic', 'kernel': 'trmm_utcopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'dir': 'generic', 'kernel': 'trmm_utcopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'dir': 'generic', 'kernel': 'trmm_ltcopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'dir': 'generic', 'kernel': 'trmm_ltcopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
        }},
    },
  },
  { 'base': '?trsm',
    'modes': {
        's': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'dir': 'generic', 'kernel': 'trsm_uncopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'dir': 'generic', 'kernel': 'trsm_uncopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'dir': 'generic', 'kernel': 'trsm_lncopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'dir': 'generic', 'kernel': 'trsm_lncopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'dir': 'generic', 'kernel': 'trsm_utcopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'dir': 'generic', 'kernel': 'trsm_utcopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'dir': 'generic', 'kernel': 'trsm_ltcopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'dir': 'generic', 'kernel': 'trsm_ltcopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'dir': 'generic', 'kernel': 'trsm_uncopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'dir': 'generic', 'kernel': 'trsm_uncopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'dir': 'generic', 'kernel': 'trsm_lncopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'dir': 'generic', 'kernel': 'trsm_lncopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'dir': 'generic', 'kernel': 'trsm_utcopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'dir': 'generic', 'kernel': 'trsm_utcopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'dir': 'generic', 'kernel': 'trsm_ltcopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'dir': 'generic', 'kernel': 'trsm_ltcopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
        }},
    },
  },
  { 'base': '?symm',
    'modes': {
        's': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iutcopy': {'dir': 'generic', 'kernel': 'symm_ucopy_8.c',
                          'addl': ['-UOUTER', '-ULOWER']},
            '_iltcopy': {'dir': 'generic', 'kernel': 'symm_lcopy_8.c',
                          'addl': ['-UOUTER', '-DLOWER']},
            '_outcopy': {'dir': 'generic', 'kernel': 'symm_ucopy_4.c',
                          'addl': ['-DOUTER', '-ULOWER']},
            '_oltcopy': {'dir': 'generic', 'kernel': 'symm_lcopy_4.c',
                          'addl': ['-DOUTER', '-DLOWER']},
        }},
    },
  },
  { 'base': '?gemm_small_kernel',
    'modes': {
        's': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tt.c',
                       },
            }
        },
        'd': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tt.c',
                       },
            }
        },
        'c': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_nc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_tc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_ct': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_cc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
            }
        },
        'z': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_nc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_tc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_ct': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_cc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
            }
        },
    },
  },
  { 'base': '?gemm_small_kernel_b0',
    'modes': {
        's': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tt.c',
                       },
            }
        },
        'd': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'gemm_small_matrix_kernel_tt.c',
                       },
            }
        },
        'c': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_nc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_tc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_ct': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_cc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
            }
        },
        'z': {
            'exts': {
                '_nn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_rr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nn.c',
                       },
                '_nt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_nc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_rc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_nt.c',
                       },
                '_tn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cn': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_cr': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tn.c',
                       },
                '_tt': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_tc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_ct': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
                '_cc': {'dir': 'generic',
                        'kernel': 'zgemm_small_matrix_kernel_tt.c',
                       },
            }
        },
    },
  },
]

kernel_confs = []
foreach _kop : base_kops
  base = _kop['base']
  modes = _kop['modes']
  # Generate the symbol flags
  _ckop_args = []
  if symb_defs.has_key(base)
    symb_base = symb_defs[base]
    if symb_base.has_key('def')
      foreach _d : symb_base['def']
        _ckop_args += ('-D' + _d)
      endforeach
    endif
    if symb_base.has_key('undef')
      foreach _u : symb_base['undef']
        _ckop_args += ('-U' + _u)
      endforeach
    endif
  endif
  foreach mode, details : modes
    # Generally, one list is required for each foreach
    __cargs = _cargs + _ckop_args
    prec_mode = precision_mappings[mode]
    # Generate the mapping for the type
    if prec_mode.has_key('def')
      foreach _d : prec_mode['def']
        __cargs += ('-D' + _d)
      endforeach
    endif
    if prec_mode.has_key('undef')
      foreach _u : prec_mode['undef']
        __cargs += ('-U' + _u)
      endforeach
    endif
    # Now the rest, one run for each ext, to get the final symbols
    foreach ext, extdat : details['exts']
      _ext_cargs = [] # Will be wiped for each ext preventing redefinitions
       # Check ext_mappings first
        if ext_mappings.has_key(ext) and (not ext_mappings.has_key('except') or base not in ext_mappings['except'])
          extmap = ext_mappings[ext]
          if extmap.has_key('def')
            foreach _d : extmap['def']
              _ext_cargs += ['-D' + _d]
            endforeach
          endif
          if extmap.has_key('undef')
            foreach _u : extmap['undef']
              _ext_cargs += ['-U' + _u]
            endforeach
          endif
        else
          # Fallback to ext_mappings_l2
          foreach ext_map : ext_mappings_l2 + ext_mappings_l3
            if ext_map['ext'] == ext and mode in ext_map['for'] and (not ext_map.has_key('except') or base not in ext_map['except'])
              if ext_map.has_key('def')
                foreach _d : ext_map['def']
                  _ext_cargs += ['-D' + _d]
                endforeach
              endif
              if ext_map.has_key('undef')
                foreach _u : ext_map['undef']
                  _ext_cargs += ['-U' + _u]
                endforeach
              endif
              break
            endif
          endforeach
        endif

      src = join_paths(extdat['dir'], extdat['kernel'])
      if extdat.has_key('addl')
        _ext_cargs += extdat['addl']
      endif
      sym_name = base.replace('?', mode) + ext
      sym_underscored = f'@sym_name@_'
      _ext_cargs += [
        f'-DASMNAME=@sym_name@',
        f'-DASMFNAME=@sym_underscored@',
        f'-DNAME=@sym_underscored@',
        f'-DCNAME=@sym_name@',
        f'-DCHAR_NAME="@sym_underscored@"',
        f'-DCHAR_CNAME="@sym_name@"',
      ]
      current_def = {
        'c_args': __cargs + _ext_cargs,
        'name': sym_name,
        'src': src
      }
      kernel_confs += current_def
    endforeach
  endforeach
endforeach

_kern_libs = []
foreach conf: kernel_confs
    message(conf['name'])
    message(conf)
    _kern_libs += static_library(
        conf['name'],
        conf['src'],
        include_directories:  _inc,
        c_args: conf['c_args'],
    )
endforeach

_kern = static_library('_kern',
                       link_whole: _kern_libs)
