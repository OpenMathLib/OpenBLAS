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
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'srot.c', 'addl': [fma3_flag]}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'drot.c', 'addl': [fma3_flag]}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zrot.c'}}},
        'cs': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zrot_sse.S'}}},
        'zd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zrot_sse2.S'}}},
        # 'xq': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zrot.S'}}},
    },
  },
  { 'base': '?swap',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'swap_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'swap_sse2.S'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zswap_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zswap_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'swap.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zswap.S'}}},
    },
  },
  { 'base': '?scal',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'sscal.c'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'dscal.c'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'cscal.c'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zscal.c'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'scal.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zscal.S'}}},
    },
  },
  { 'base': '?copy',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'copy_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'copy_sse2.S'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zcopy_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zcopy_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'copy.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zcopy.S'}}},
    },
  },
  { 'base': '?axpy',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'saxpy.c'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'daxpy.c'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'caxpy.c'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zaxpy.c'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'axpy.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zaxpy.S'}}},
    },
  },
  { 'base': '?dot',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'sdot.c'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'ddot.c'}}},
    },
  },
  { 'base': '?dotc',
    'modes': {
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'cdot.c'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zdot.c'}}},
    },
  },
  { 'base': '?dotu',
    'modes': {
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'cdot.c'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zdot.c'}}},
    },
  },
  # TODO(rg): Check?
  { 'base': '?dsdot',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'sdot.c'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'sdot.c'}}},
    },
  },
  # TODO(rg): Add dsdotkernel conditionals
  # xDOTU xDOTC xxDOT aren't present
  { 'base': '?nrm2',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'nrm2_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'nrm2.S'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'znrm2_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'znrm2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'nrm2.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'znrm2.S'}}},
    },
  },
  { 'base': '?asum',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'sasum.c'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'dasum.c'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zasum_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zasum_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'asum.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zasum.S'}}},
    },
  },
  { 'base': '?amax',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse2.S'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zamax_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zamax_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zamax.S'}}},
    },
  },
  { 'base': '?sum',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'sum.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'sum.c'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zsum_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zsum_sse2.S'}}},
    },
  },
  { 'base': '?amin',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse2.S'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zamax_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'zamax_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'izamax.S'}}},
    },
  },
  { 'base': 'i?amax',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax_sse2.S'}}},
        'c': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'izamax_sse.S'}}},
        'z': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'izamax_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
        # 'x': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'izamax.S'}}},
    },
  },
  { 'base': 'i?amin',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'iamin.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'izamin.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'izamin.c'}}},
    },
  },
  { 'base': 'i?max',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
    },
  },
  { 'base': 'i?min',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'iamax.S'}}},
    },
  },
  { 'base': '?max',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax.S'}}},
    },
  },
  { 'base': '?min',
    'modes': {
        's': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse.S'}}},
        'd': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax_sse2.S'}}},
        # 'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'amax.S'}}},
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
                '_n': {'dir': 'x86_64', 'kernel': 'sgemv_n_4.c'},
                '_t': {'dir': 'x86_64', 'kernel': 'sgemv_t_4.c'},
            }
        },
        'd': {
            'exts': {
                '_n': {'dir': 'x86_64', 'kernel': 'dgemv_n_4.c'},
                '_t': {'dir': 'x86_64', 'kernel': 'dgemv_t_4.c'},
            }
        },
        # 'q': {
        #     'exts': {
        #         '_n': {'dir': 'x86_64', 'kernel': 'qgemv_n.S'},
        #         '_t': {'dir': 'x86_64', 'kernel': 'qgemv_t.S'},
        #     }
        # },
        'c': {
            'exts': {
                '_n': {'dir': 'x86_64', 'kernel': 'cgemv_n_4.c'},
                '_t': {'dir': 'x86_64', 'kernel': 'cgemv_t_4.c'},
                '_r': {'dir': 'x86_64', 'kernel': 'cgemv_n_4.c'},
                '_c': {'dir': 'x86_64', 'kernel': 'cgemv_t_4.c'},
                '_o': {'dir': 'x86_64', 'kernel': 'cgemv_n_4.c'},
                '_u': {'dir': 'x86_64', 'kernel': 'cgemv_t_4.c'},
                '_s': {'dir': 'x86_64', 'kernel': 'cgemv_n_4.c'},
                '_d': {'dir': 'x86_64', 'kernel': 'cgemv_t_4.c'},
            }
        },
        'z': {
            'exts': {
                '_n': {'dir': 'x86_64', 'kernel': 'zgemv_n_4.c'},
                '_t': {'dir': 'x86_64', 'kernel': 'zgemv_t_4.c'},
                '_r': {'dir': 'x86_64', 'kernel': 'zgemv_n_4.c'},
                '_c': {'dir': 'x86_64', 'kernel': 'zgemv_t_4.c'},
                '_o': {'dir': 'x86_64', 'kernel': 'zgemv_n_4.c'},
                '_u': {'dir': 'x86_64', 'kernel': 'zgemv_t_4.c'},
                '_s': {'dir': 'x86_64', 'kernel': 'zgemv_n_4.c'},
                '_d': {'dir': 'x86_64', 'kernel': 'zgemv_t_4.c'},
            }
        },
        # 'x': {
        #     'exts': {
        #         '_n': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #         '_t': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #         '_r': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #         '_c': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #         '_o': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #         '_u': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #         '_s': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #         '_d': {'dir': 'arm', 'kernel': 'zgemv_n.c'},
        #     }
        # },
    },
  },
  { 'base': '?symv',
    'modes': {
        's': {
            'exts': {
                '_U': {'dir': 'x86_64', 'kernel': 'ssymv_U.c'},
                '_L': {'dir': 'x86_64', 'kernel': 'ssymv_L.c'},
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
                '_U': {'dir': 'x86_64', 'kernel': 'zsymv_U_sse2.S'},
                '_L': {'dir': 'x86_64', 'kernel': 'zsymv_L_sse2.S'},
            }
        },
        # 'q': {
        #     'exts': {
        #         '_U': {'dir': 'generic', 'kernel': 'symv_k.c'},
        #         '_L': {'dir': 'generic', 'kernel': 'symv_k.c'},
        #     }
        # },
        # 'x': {
        #     'exts': {
        #         '_U': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
        #         '_L': {'dir': 'generic', 'kernel': 'zsymv_k.c'},
        #     }
        # },
    },
  },
  { 'base': '?lsame',
    'modes': {
        '': {'exts': {'': {'dir': 'x86_64', 'kernel': 'lsame.S'}}},
    },
  },
  { 'base': '?cabs',
    'modes': {
        's': {'exts': {'1': {'dir': 'x86_64', 'kernel': 'cabs.S'}}},
        'd': {'exts': {'1': {'dir': 'x86_64', 'kernel': 'cabs.S'}}},
        # 'q': {'exts': {'': {'dir': 'generic', 'kernel': 'cabs.c'}}},
    },
  },
  { 'base': '?gemm3m',
    'modes': {
        'c': {'exts': {
            '_kernel': {'dir': 'x86_64', 'kernel': 'cgemm3m_kernel_8x4_haswell.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_oncopyb': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_4.c',
                        'addl': ['-DUSE_ALPHA']},
            '_otcopyb': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_4.c',
                        'addl': ['-DUSE_ALPHA']},
            '_itcopyb': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_itcopyr': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_itcopyi': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_incopyb': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_incopyr': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_oncopyr': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_otcopyr': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_incopyi': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_oncopyi': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DIMAGE_ONLY']},
            '_otcopyi': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
        }},
        'z': {'exts': {
            '_kernel': {'dir': 'x86_64', 'kernel': 'zgemm3m_kernel_4x4_haswell.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_oncopyb': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_4.c',
                        'addl': ['-DUSE_ALPHA']},
            '_otcopyb': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_4.c',
                        'addl': ['-DUSE_ALPHA']},
            '_itcopyb': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_itcopyr': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_itcopyi': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_incopyb': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_incopyr': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_oncopyr': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_otcopyr': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_incopyi': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_8.c',
                        'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_oncopyi': {'dir': 'generic', 'kernel': 'zgemm3m_ncopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DIMAGE_ONLY']},
            '_otcopyi': {'dir': 'generic', 'kernel': 'zgemm3m_tcopy_4.c',
                        'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
        }},
    },
  },
  { 'base': '?ger',
    'modes': {
        's': {'exts': {'_k': {'dir': 'generic', 'kernel': 'ger.c'}}},
        'd': {'exts': {'_k': {'dir': 'generic', 'kernel': 'ger.c'}}},
        # 'q': {'exts': {'_k': {'dir': 'generic', 'kernel': 'ger.c'}}},
    },
  },
  { 'base': '?geru',
    'modes': {
        'c': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        'z': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        # 'x': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
    },
  },
  { 'base': '?gerc',
    'modes': {
        'c': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        'z': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
        # 'x': {'exts': {'_k': {'dir': 'generic', 'kernel': 'zger.c'}}},
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
                '_U': {'dir': 'x86_64', 'kernel': 'zsymv_U_sse2.S'},
                '_L': {'dir': 'x86_64', 'kernel': 'zsymv_L_sse2.S'},
                '_V': {'dir': 'x86_64', 'kernel': 'zsymv_L_sse2.S'},
                '_M': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
            }
        },
        # 'x': {
        #     'exts': {
        #         '_U': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
        #         '_L': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
        #         '_V': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
        #         '_M': {'dir': 'generic', 'kernel': 'zhemv_k.c'},
        #     }
        # },
    },
  },
  # { 'base': '?bgemv',
  #   'modes': {
  #       's': {
  #           'exts': {
  #               '_n': {'dir': 'x86_64', 'kernel': 'sbgemv_n.c'},
  #               '_t': {'dir': 'x86_64', 'kernel': 'sbgemv_n.c'},
  #           }
  #       }
  #   },
  # },
  # Level 3 symbols
  { 'base': '?gemm_kernel',
    'modes': {
        's': {'exts': {'': {'dir': 'x86_64', 'kernel': 'sgemm_kernel_8x4_haswell_2.c'}}},
        'd': {'exts': {'': {'dir': 'x86_64', 'kernel': 'dgemm_kernel_4x8_haswell.S'}}},
        'c': {
            'exts': {
                '_n': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.c', 'addl': ['-DNN']},
                '_l': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.c', 'addl': ['-DCN']},
                # TODO(rg): What about _r conditionals? Makefile.L3:2969
                '_r': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.c', 'addl': ['-DNC']},
                '_b': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.c', 'addl': ['-DCC']},
            }
        },
        'z': {
            'exts': {
                '_n': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.c', 'addl': ['-DNN']},
                '_l': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.c', 'addl': ['-DCN']},
                '_r': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.c', 'addl': ['-DNC']},
                '_b': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.c', 'addl': ['-DCC']},
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
                '_LN': {'dir': 'x86_64', 'kernel': 'sgemm_kernel_8x4_haswell.c'},
                '_LT': {'dir': 'x86_64', 'kernel': 'sgemm_kernel_8x4_haswell.c', 'addl': ['-DLEFT', '-DTRANSA']},
                '_RN': {'dir': 'x86_64', 'kernel': 'sgemm_kernel_8x4_haswell.c'},
                '_RT': {'dir': 'x86_64', 'kernel': 'sgemm_kernel_8x4_haswell.c'},
            }
        },
        'd': {
            'exts': {
                '_LN': {'dir': 'x86_64', 'kernel': 'dtrmm_kernel_4x8_haswell.c'},
                '_LT': {'dir': 'x86_64', 'kernel': 'dtrmm_kernel_4x8_haswell.c', 'addl': ['-DLEFT', '-DTRANSA']},
                '_RN': {'dir': 'x86_64', 'kernel': 'dtrmm_kernel_4x8_haswell.c'},
                '_RT': {'dir': 'x86_64', 'kernel': 'dtrmm_kernel_4x8_haswell.c'},
            }
        },
        'c': {
            'exts': {
                '_LN': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-UCONJ', '-DNN']},
                '_LT': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-DLEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_LR': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-DCONJ', '-DCN']},
                '_LC': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-DCONJ', '-DCN']},
                '_RN': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-UCONJ', '-DNN']},
                '_RT': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-ULEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_RR': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-DCONJ', '-DNC']},
                '_RC': {'dir': 'x86_64', 'kernel': 'cgemm_kernel_8x2_haswell.S',
                        'addl': ['-DCONJ', '-DCN']},
            }
        },
        'z': {
            'exts': {
                '_LN': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-UCONJ', '-DNN']},
                '_LT': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-DLEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_LR': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-DCONJ', '-DCN']},
                '_LC': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-DCONJ', '-DCN']},
                '_RN': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-UCONJ', '-DNN']},
                '_RT': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-ULEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_RR': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-DCONJ', '-DNC']},
                '_RC': {'dir': 'x86_64', 'kernel': 'zgemm_kernel_4x2_haswell.S',
                        'addl': ['-DCONJ', '-DCN']},
            },
        },
    },
  },
  { 'base': '?trsm_kernel',
    'modes': {
        's': {
            'exts': {
                '_LN': {'dir': 'x86_64', 'kernel': 'strsm_kernel_8x4_haswell_LN.c',
                        'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'dir': 'x86_64', 'kernel': 'strsm_kernel_8x4_haswell_LT.c',
                        'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_RN': {'dir': 'x86_64', 'kernel': 'strsm_kernel_8x4_haswell_RN.c',
                        'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'dir': 'x86_64', 'kernel': 'strsm_kernel_8x4_haswell_RT.c',
                        'addl': ['-DRT', '-UUPPER', '-UCONJ']},
            },
        },
        'd': {
            'exts': {
                '_LN': {'dir': 'generic', 'kernel': 'trsm_kernel_LN.c',
                        'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'dir': 'generic', 'kernel': 'trsm_kernel_LT.c',
                        'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_RN': {'dir': 'x86_64', 'kernel': 'dtrsm_kernel_RN_haswell.c',
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
        'z': {
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
  { 'base': '?gemm',
    'modes': {
        's': {'exts': {
            '_beta': {'dir': 'x86_64', 'kernel': 'sgemm_beta_skylakex.c'},
            '_small_matrix_permit': {'dir': 'generic',
                                     'kernel': 'gemm_small_matrix_permit.c'},
            # TODO(rg): the _NUM prefixes are arch dependent
            '_incopy': {'dir': 'generic', 'kernel': 'gemm_ncopy_8.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_8.c'},
            '_oncopy': {'dir': 'x86_64', 'kernel': 'sgemm_ncopy_4_skylakex.c'},
            '_otcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_4.c'},
            # TODO(rg): direct and direct_performant are built only conditionally
            '_direct': {'dir': 'x86_64', 'kernel': 'sgemm_direct_skylakex.c'},
            '_direct_performant': {'dir': 'x86_64', 'kernel': 'sgemm_direct_performant.c'},
        }},
        'd': {'exts': {
            '_beta': {'dir': 'x86_64', 'kernel': 'dgemm_beta_skylakex.c'},
            '_small_matrix_permit': {'dir': 'generic',
                                     'kernel': 'gemm_small_matrix_permit.c'},
            '_incopy': {'dir': 'generic', 'kernel': 'gemm_ncopy_4.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_4.c'},
            '_oncopy': {'dir': 'x86_64', 'kernel': 'dgemm_ncopy_8_skylakex.c'},
            '_otcopy': {'dir': 'generic', 'kernel': 'gemm_tcopy_8.c'},
        }},
        'c': {'exts': {
            '_beta': {'dir': 'x86_64', 'kernel': 'zgemm_beta.S'},
            '_small_matrix_permit': {'dir': 'generic',
                                     'kernel': 'zgemm_small_matrix_permit.c'},
            '_incopy': {'dir': 'generic', 'kernel': 'zgemm_ncopy_8.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'zgemm_tcopy_8.c'},
            '_oncopy': {'dir': 'generic', 'kernel': 'zgemm_ncopy_2.c'},
            '_otcopy': {'dir': 'generic', 'kernel': 'zgemm_tcopy_2.c'},
        }},
        'z': {'exts': {
            '_beta': {'dir': 'x86_64', 'kernel': 'zgemm_beta.S'},
            '_small_matrix_permit': {'dir': 'generic',
                                     'kernel': 'zgemm_small_matrix_permit.c'},
            '_incopy': {'dir': 'generic', 'kernel': 'zgemm_ncopy_4.c'},
            '_itcopy': {'dir': 'generic', 'kernel': 'zgemm_tcopy_4.c'},
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
        'd': {'exts': {
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
        'q': {'exts': {
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
        'c': {'exts': {
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
        'z': {'exts': {
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
  { 'base': '?omatcopy_k',
    'modes': {
        's': {'exts': {
            '_cn': {'dir': 'arm', 'kernel': 'omatcopy_cn.c',
                          'addl': ['-UROWM']},
            '_rn': {'dir': 'arm', 'kernel': 'omatcopy_rn.c',
                          'addl': ['-DROWM']},
            '_ct': {'dir': 'arm', 'kernel': 'omatcopy_ct.c',
                          'addl': ['-UROWM']},
            '_rt': {'dir': 'x86_64', 'kernel': 'omatcopy_rt.c',
                          'addl': ['-DROWM']},
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
                # '_b0_nn': {'dir': 'generic',
                #            'kernel': 'gemm_small_matrix_kernel_nn.c',
                #            'addl': ['-DB0'],
                #        },
                # '_b0_nt': {'dir': 'generic',
                #            'kernel': 'gemm_small_matrix_kernel_nt.c',
                #            'addl': ['-DB0'],
                #        },
                # '_b0_tn': {'dir': 'generic',
                #            'kernel': 'gemm_small_matrix_kernel_tn.c',
                #            'addl': ['-DB0'],
                #        },
                # '_b0_tt': {'dir': 'generic',
                #            'kernel': 'gemm_small_matrix_kernel_tt.c',
                #            'addl': ['-DB0'],
                #        },
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
    if mode == 'x' or mode == 'q'
      continue
    endif
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
_kern_deps = []
_is_asm = false
foreach conf: kernel_confs
  _kern_libs += static_library(
      conf['name'],
      conf['src'],
      include_directories:  _inc,
      c_args: conf['c_args'],
      # See gh discussion 13374 for why, basically .S are coded as fortran..
      fortran_args: conf['c_args'],
  )
endforeach

_kern = static_library('_kern',
                       link_whole: _kern_libs,
                       dependencies: _kern_deps)
