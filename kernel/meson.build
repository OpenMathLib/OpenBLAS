# Ordered As per https://netlib.org/blas/blasqr.pdf
# NOTE: xROTG xROTMG xROTM have no kernels?
# TODO: Actually test and set this
fma3_flag = []
if conf_hdat.has('HAVE_FMA3')
  fma3_flag += '-mfma'
endif

SGEMM_UNROLL_N = makefile_conf_dat.get('SGEMM_UNROLL_N')
SGEMM_UNROLL_M = makefile_conf_dat.get('SGEMM_UNROLL_M')
DGEMM_UNROLL_N = makefile_conf_dat.get('DGEMM_UNROLL_N')
DGEMM_UNROLL_M = makefile_conf_dat.get('DGEMM_UNROLL_M')
CGEMM_UNROLL_N = makefile_conf_dat.get('CGEMM_UNROLL_N')
CGEMM_UNROLL_M = makefile_conf_dat.get('CGEMM_UNROLL_M')
ZGEMM_UNROLL_N = makefile_conf_dat.get('ZGEMM_UNROLL_N')
ZGEMM_UNROLL_M = makefile_conf_dat.get('ZGEMM_UNROLL_M')
CGEMM3M_UNROLL_N = makefile_conf_dat.get('CGEMM3M_UNROLL_N')
CGEMM3M_UNROLL_M = makefile_conf_dat.get('CGEMM3M_UNROLL_M')
ZGEMM3M_UNROLL_N = makefile_conf_dat.get('ZGEMM3M_UNROLL_N')
ZGEMM3M_UNROLL_M = makefile_conf_dat.get('ZGEMM3M_UNROLL_M')

subdir('meson_base')
subdir('x86_64')
subdir('arm64')

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
        's': {'exts': {'_k': {'addl': fma3_flag}}},
        'd': {'exts': {'_k': {'addl': fma3_flag}}},
        # 'q': {'exts': {'_k': {}}},
        'cs': {'exts': {'_k': {}}},
        'zd': {'exts': {'_k': {}}},
        # 'xq': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?swap',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?scal',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?copy',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?axpy',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {'addl': ['-UCONJ']}}},
        'z': {'exts': {'_k': {'addl': ['-UCONJ']}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?axpyc',
    'modes': {
        'c': {'exts': {'_k': {'addl': ['-DCONJ']}}},
        'z': {'exts': {'_k': {'addl': ['-DCONJ']}}},
    },
  },
  { 'base': '?dot',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?dotc',
    'modes': {
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?dotu',
    'modes': {
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
    },
  },
  # TODO(rg): Check?
  { 'base': '?dsdot',
    'modes': {
        's': {'exts': {'_k': {}}},
        '': {'exts': {'_k': {'addl': ['-DDSDOT']}}},
    },
  },
  # TODO(rg): Add dsdotkernel conditionals
  # xDOTU xDOTC xxDOT aren't present
  { 'base': '?nrm2',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?asum',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?amax',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?sum',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?amin',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': 'i?amax',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': 'i?amin',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
    },
  },
  { 'base': 'i?max',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
    },
  },
  { 'base': 'i?min',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?max',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?min',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?axpby',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        'c': {'exts': {'_k': {}}},
        'z': {'exts': {'_k': {}}},
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
                '_n': {},
                '_t': {},
            }
        },
        'd': {
            'exts': {
                '_n': {},
                '_t': {},
            }
        },
        # 'q': {
        #     'exts': {
        #         '_n': {},
        #         '_t': {},
        #     }
        # },
        'c': {
            'exts': {
                '_n': {},
                '_t': {},
                '_r': {},
                '_c': {},
                '_o': {},
                '_u': {},
                '_s': {},
                '_d': {},
            }
        },
        'z': {
            'exts': {
                '_n': {},
                '_t': {},
                '_r': {},
                '_c': {},
                '_o': {},
                '_u': {},
                '_s': {},
                '_d': {},
            }
        },
        # 'x': {
        #     'exts': {
        #         '_n': {},
        #         '_t': {},
        #         '_r': {},
        #         '_c': {},
        #         '_o': {},
        #         '_u': {},
        #         '_s': {},
        #         '_d': {},
        #     }
        # },
    },
  },
  { 'base': '?symv',
    'modes': {
        's': {
            'exts': {
                '_U': {},
                '_L': {},
            }
        },
        'd': {
            'exts': {
                '_U': {},
                '_L': {},
            }
        },
        'c': {
            'exts': {
                '_U': {},
                '_L': {},
            }
        },
        'z': {
            'exts': {
                '_U': {},
                '_L': {},
            }
        },
        # 'q': {
        #     'exts': {
        #         '_U': {},
        #         '_L': {},
        #     }
        # },
        # 'x': {
        #     'exts': {
        #         '_U': {},
        #         '_L': {},
        #     }
        # },
    },
  },
  { 'base': '?lsame',
    'modes': {
        '': {'exts': {'': {}}},
    },
  },
  { 'base': '?cabs',
    'modes': {
        's': {'exts': {'1': {}}},
        'd': {'exts': {'1': {}}},
        # 'q': {'exts': {'': {}}},
    },
  },
  { 'base': '?gemm3m',
    'modes': {
        'c': {'exts': {
            '_kernel': {'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_oncopyb': {'addl': ['-DUSE_ALPHA']},
            '_otcopyb': {'addl': ['-DUSE_ALPHA']},
            '_itcopyb': {'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_itcopyr': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_itcopyi': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_incopyb': {'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_incopyr': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_oncopyr': {'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_otcopyr': {'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_incopyi': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_oncopyi': {'addl': ['-DUSE_ALPHA', '-DIMAGE_ONLY']},
            '_otcopyi': {'addl': ['-DUSE_ALPHA', '-DIMAGE_ONLY']},
        }},
        'z': {'exts': {
            '_kernel': {'addl': ['-DNN']},
            '_oncopyb': {'addl': ['-DUSE_ALPHA']},
            '_otcopyb': {'addl': ['-DUSE_ALPHA']},
            '_itcopyb': {'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_itcopyr': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_itcopyi': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_incopyb': {'addl': ['-DICOPY', '-UUSE_ALPHA']},
            '_incopyr': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DREAL_ONLY']},
            '_oncopyr': {'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_otcopyr': {'addl': ['-DUSE_ALPHA', '-DREAL_ONLY']},
            '_incopyi': {'addl': ['-DICOPY', '-UUSE_ALPHA', '-DIMAGE_ONLY']},
            '_oncopyi': {'addl': ['-DUSE_ALPHA', '-DIMAGE_ONLY']},
            '_otcopyi': {'addl': ['-DUSE_ALPHA', '-DIMAGE_ONLY']},
        }},
    },
  },
  { 'base': '?ger',
    'modes': {
        's': {'exts': {'_k': {}}},
        'd': {'exts': {'_k': {}}},
        # 'q': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?geru',
    'modes': {
        'c': {'exts': {'_k': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UCONJ']}}},
        'z': {'exts': {'_k': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UCONJ']}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?gerc',
    'modes': {
        'c': {'exts': {'_k': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DCONJ']}}},
        'z': {'exts': {'_k': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DCONJ']}}},
        # 'x': {'exts': {'_k': {}}},
    },
  },
  { 'base': '?gerv',
    'modes': {
        'c': {'exts': {'_k': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UCONJ', '-DXCONJ']}}},
        'z': {'exts': {'_k': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UCONJ', '-DXCONJ']}}},
    },
  },
  { 'base': '?hemv',
    'modes': {
        'c': {
            'exts': {
                '_U': {'addl': ['-DHEMV']},
                '_L': {'addl': ['-DHEMV']},
                '_V': {'addl': ['-DHEMV', '-DHEMVREV']},
                '_M': {'addl': ['-DHEMV', '-DHEMVREV']},
            }
        },
        'z': {
            'exts': {
                '_U': {'addl': ['-DHEMV']},
                '_L': {'addl': ['-DHEMV']},
                '_V': {'addl': ['-DHEMV', '-DHEMVREV']},
                '_M': {'addl': ['-DHEMV', '-DHEMVREV']},
            }
        },
        # 'x': {
        #     'exts': {
        #         '_U': {},
        #         '_L': {},
        #         '_V': {},
        #         '_M': {},
        #     }
        # },
    },
  },
  # { 'base': '?bgemv',
  #   'modes': {
  #       's': {
  #           'exts': {
  #               '_n': {},
  #               '_t': {},
  #           }
  #       }
  #   },
  # },
  # Level 3 symbols
  { 'base': '?gemm_kernel',
    'modes': {
        's': {'exts': {'': {}}},
        'd': {'exts': {'': {}}},
        'c': {
            'exts': {
                '_n': {'addl': ['-DNN']},
                '_l': {'addl': ['-DCN']},
                # TODO(rg): What about _r conditionals? Makefile.L3:2969
                '_r': {'addl': ['-DNC']},
                '_b': {'addl': ['-DCC']},
            }
        },
        'z': {
            'exts': {
                '_n': {'addl': ['-DNN']},
                '_l': {'addl': ['-DCN']},
                '_r': {'addl': ['-DNC']},
                '_b': {'addl': ['-DCC']},
            }
        }
        # 'q': {'exts': {'': {}}},
        # 'x': {'exts': {'': {}}},
    },
  },
  { 'base': '?trmm_kernel',
    'modes': {
        's': {
            'exts': {
                '_LN': {},
                '_LT': {'addl': ['-DLEFT', '-DTRANSA']},
                '_RN': {},
                '_RT': {},
            }
        },
        'd': {
            'exts': {
                '_LN': {},
                '_LT': {'addl': ['-DLEFT', '-DTRANSA']},
                '_RN': {},
                '_RT': {},
            }
        },
        'c': {
            'exts': {
                '_LN': {'addl': ['-DLEFT', '-UTRANSA', '-UCONJ', '-DNN']},
                '_LT': {'addl': ['-DLEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_LR': {'addl': ['-DLEFT', '-UTRANSA', '-DCONJ', '-DCN']},
                '_LC': {'addl': ['-DLEFT', '-DTRANSA', '-DCONJ', '-DCN']},
                '_RN': {'addl': ['-ULEFT', '-UTRANSA', '-UCONJ', '-DNN']},
                '_RT': {'addl': ['-ULEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_RR': {'addl': ['-ULEFT', '-UTRANSA', '-DCONJ', '-DNC']},
                '_RC': {'addl': ['-ULEFT', '-DTRANSA', '-DCONJ', '-DNC']},
            }
        },
        'z': {
            'exts': {
                '_LN': {'addl': ['-DLEFT', '-UTRANSA', '-UCONJ', '-DNN']},
                '_LT': {'addl': ['-DLEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_LR': {'addl': ['-DLEFT', '-UTRANSA', '-DCONJ', '-DCN']},
                '_LC': {'addl': ['-DLEFT', '-DTRANSA', '-DCONJ', '-DCN']},
                '_RN': {'addl': ['-ULEFT', '-UTRANSA', '-UCONJ', '-DNN']},
                '_RT': {'addl': ['-ULEFT', '-DTRANSA', '-UCONJ', '-DNN']},
                '_RR': {'addl': ['-ULEFT', '-UTRANSA', '-DCONJ', '-DNC']},
                '_RC': {'addl': ['-ULEFT', '-DTRANSA', '-DCONJ', '-DNC']},
            },
        },
    },
  },
  { 'base': '?trsm_kernel',
    'modes': {
        's': {
            'exts': {
                '_LN': {'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_RN': {'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'addl': ['-DRT', '-UUPPER', '-UCONJ']},
            },
        },
        'd': {
            'exts': {
                '_LN': {'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_RN': {'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'addl': ['-DRT', '-UUPPER', '-UCONJ']},
            },
        },
        'c': {
            'exts': {
                '_LN': {'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_LR': {'addl': ['-DLN', '-DUPPER', '-DCONJ']},
                '_LC': {'addl': ['-DLT', '-UUPPER', '-DCONJ']},
                '_RN': {'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'addl': ['-DRT', '-UUPPER', '-UCONJ']},
                '_RR': {'addl': ['-DRN', '-DUPPER', '-DCONJ']},
                '_RC': {'addl': ['-DRT', '-UUPPER', '-DCONJ']},
            },
        },
        'z': {
            'exts': {
                '_LN': {'addl': ['-DLN', '-DUPPER', '-UCONJ']},
                '_LT': {'addl': ['-DLT', '-UUPPER', '-UCONJ']},
                '_LR': {'addl': ['-DLN', '-DUPPER', '-DCONJ']},
                '_LC': {'addl': ['-DLT', '-UUPPER', '-DCONJ']},
                '_RN': {'addl': ['-DRN', '-DUPPER', '-UCONJ']},
                '_RT': {'addl': ['-DRT', '-UUPPER', '-UCONJ']},
                '_RR': {'addl': ['-DRN', '-DUPPER', '-DCONJ']},
                '_RC': {'addl': ['-DRT', '-UUPPER', '-DCONJ']},
            },
        },
    },
  },
  { 'base': '?gemm',
    'modes': {
        's': {'exts': {
            '_beta': {},
            '_small_matrix_permit': {},
            # TODO(rg): the _NUM prefixes are arch dependent
            '_incopy': {},
            '_itcopy': {},
            '_oncopy': {},
            '_otcopy': {},
            # TODO(rg): direct and direct_performant are built only conditionally
            '_direct': {},
            '_direct_performant': {},
        }},
        'd': {'exts': {
            '_beta': {},
            '_small_matrix_permit': {},
            '_incopy': {},
            '_itcopy': {},
            '_oncopy': {},
            '_otcopy': {},
        }},
        'c': {'exts': {
            '_beta': {},
            '_small_matrix_permit': {},
            '_incopy': {},
            '_itcopy': {},
            '_oncopy': {},
            '_otcopy': {},
        }},
        'z': {'exts': {
            '_beta': {},
            '_small_matrix_permit': {},
            '_incopy': {},
            '_itcopy': {},
            '_oncopy': {},
            '_otcopy': {},
        }},
    },
  },
  { 'base': '?trmm',
    'modes': {
        's': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
        }},
        'd': {'exts': {
            '_iunucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
        }},
        'c': {'exts': {
            '_iunucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
        }},
        'z': {'exts': {
            '_iunucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER', '-UUNIT']},
        }},
    },
  },
  { 'base': '?hemm',
    'modes': {
        'c': {'exts': {
            '_iutcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER']},
            '_iltcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER']},
            '_outcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER']},
            '_oltcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER']},
        }},
        'z': {'exts': {
            '_iutcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER']},
            '_iltcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER']},
            '_outcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER']},
            '_oltcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER']},
        }},
    },
  },
  { 'base': '?trsm',
    'modes': {
        's': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
        }},
        'd': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
        }},
        'q': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
        }},
        'c': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
        }},
        'z': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iunucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iunncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_ilnucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_ilnncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_iutucopy': {'addl': ['-UOUTER', '-ULOWER', '-DUNIT']},
            '_iutncopy': {'addl': ['-UOUTER', '-ULOWER', '-UUNIT']},
            '_iltucopy': {'addl': ['-UOUTER', '-DLOWER', '-DUNIT']},
            '_iltncopy': {'addl': ['-UOUTER', '-DLOWER', '-UUNIT']},
            '_ounucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_ounncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_olnucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_olnncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
            '_outucopy': {'addl': ['-DOUTER', '-ULOWER', '-DUNIT']},
            '_outncopy': {'addl': ['-DOUTER', '-ULOWER', '-UUNIT']},
            '_oltucopy': {'addl': ['-DOUTER', '-DLOWER', '-DUNIT']},
            '_oltncopy': {'addl': ['-DOUTER', '-DLOWER', '-UUNIT']},
        }},
    },
  },
  { 'base': '?symm',
    'modes': {
        's': {'exts': {
            # TODO(rg): These actually use $(SGEMM_UNROLL_M) to choose the size
            '_iutcopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER']},
            '_iltcopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER']},
            '_outcopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER']},
            '_oltcopy': {'addl': ['-UDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER']},
        }},
        'd': {'exts': {
            '_iutcopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-ULOWER']},
            '_iltcopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-UOUTER', '-DLOWER']},
            '_outcopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-ULOWER']},
            '_oltcopy': {'addl': ['-DDOUBLE', '-UCOMPLEX', '-DOUTER', '-DLOWER']},
        }},
        'c': {'exts': {
            '_iutcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER']},
            '_iltcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER']},
            '_outcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER']},
            '_oltcopy': {'addl': ['-UDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER']},
        }},
        'z': {'exts': {
            '_iutcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-ULOWER']},
            '_iltcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-UOUTER', '-DLOWER']},
            '_outcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-ULOWER']},
            '_oltcopy': {'addl': ['-DDOUBLE', '-DCOMPLEX', '-DOUTER', '-DLOWER']},
        }},
    },
  },
  { 'base': '?omatcopy_k',
    'modes': {
        's': {'exts': {
            '_cn': {'addl': ['-UROWM']},
            '_rn': {'addl': ['-DROWM']},
            '_ct': {'addl': ['-UROWM']},
            '_rt': {'addl': ['-DROWM']},
        }},
        'd': {'exts': {
            '_cn': {'addl': ['-UROWM']},
            '_rn': {'addl': ['-DROWM']},
            '_ct': {'addl': ['-UROWM']},
            '_rt': {'addl': ['-DROWM']},
        }},
        'c': {'exts': {
            '_cn': {'addl': ['-UROWM', '-UCONJ']},
            '_rn': {'addl': ['-DROWM', '-UCONJ']},
            '_ct': {'addl': ['-UROWM', '-UCONJ']},
            '_rt': {'addl': ['-DROWM', '-UCONJ']},
            '_cnc': {'addl': ['-UROWM', '-DCONJ']},
            '_rnc': {'addl': ['-DROWM', '-DCONJ']},
            '_ctc': {'addl': ['-UROWM', '-DCONJ']},
            '_rtc': {'addl': ['-DROWM', '-DCONJ']},
        }},
        'z': {'exts': {
            '_cn': {'addl': ['-UROWM', '-UCONJ']},
            '_rn': {'addl': ['-DROWM', '-UCONJ']},
            '_ct': {'addl': ['-UROWM', '-UCONJ']},
            '_rt': {'addl': ['-DROWM', '-UCONJ']},
            '_cnc': {'addl': ['-UROWM', '-DCONJ']},
            '_rnc': {'addl': ['-DROWM', '-DCONJ']},
            '_ctc': {'addl': ['-UROWM', '-DCONJ']},
            '_rtc': {'addl': ['-DROWM', '-DCONJ']},
        }},
    },
  },
  { 'base': '?imatcopy_k',
    'modes': {
        's': {'exts': {
            '_cn': {'addl': ['-UROWM']},
            '_rn': {'addl': ['-DROWM']},
            '_ct': {'addl': ['-UROWM']},
            '_rt': {'addl': ['-DROWM']},
        }},
        'd': {'exts': {
            '_cn': {'addl': ['-UROWM']},
            '_rn': {'addl': ['-DROWM']},
            '_ct': {'addl': ['-UROWM']},
            '_rt': {'addl': ['-DROWM']},
        }},
        'c': {'exts': {
            '_cn': {'addl': ['-UROWM', '-UCONJ']},
            '_rn': {'addl': ['-DROWM', '-UCONJ']},
            '_ct': {'addl': ['-UROWM', '-UCONJ']},
            '_rt': {'addl': ['-DROWM', '-UCONJ']},
            '_cnc': {'addl': ['-UROWM', '-DCONJ']},
            '_rnc': {'addl': ['-DROWM', '-DCONJ']},
            '_ctc': {'addl': ['-UROWM', '-DCONJ']},
            '_rtc': {'addl': ['-DROWM', '-DCONJ']},
        }},
        'z': {'exts': {
            '_cn': {'addl': ['-UROWM', '-UCONJ']},
            '_rn': {'addl': ['-DROWM', '-UCONJ']},
            '_ct': {'addl': ['-UROWM', '-UCONJ']},
            '_rt': {'addl': ['-DROWM', '-UCONJ']},
            '_cnc': {'addl': ['-UROWM', '-DCONJ']},
            '_rnc': {'addl': ['-DROWM', '-DCONJ']},
            '_ctc': {'addl': ['-UROWM', '-DCONJ']},
            '_rtc': {'addl': ['-DROWM', '-DCONJ']},
        }},
    },
  },
  {
    'base': '?geadd',
    'modes': {
        's': {'exts': {'_k': {'addl': ['-UROWM']}}},
        'd': {'exts': {'_k': {'addl': ['-UROWM']}}},
        'c': {'exts': {'_k': {'addl': ['-UROWM']}}},
        'z': {'exts': {'_k': {'addl': ['-UROWM']}}},
    },
  },
  { 'base': '?gemm_small_kernel',
    'modes': {
        's': {
            'exts': {
                '_nn': {},
                '_nt': {},
                '_tn': {},
                '_tt': {},
                # '_b0_nn': {'addl': ['-DB0']},
                # '_b0_nt': {'addl': ['-DB0']},
                # '_b0_tn': {'addl': ['-DB0']},
                # '_b0_tt': {'addl': ['-DB0']},
            }
        },
        'd': {
            'exts': {
                '_nn': {},
                '_nt': {},
                '_tn': {},
                '_tt': {},
            }
        },
        'c': {
            'exts': {
                '_nn': {},
                '_nr': {},
                '_rn': {},
                '_rr': {},
                '_nt': {},
                '_nc': {},
                '_rt': {},
                '_rc': {},
                '_tn': {},
                '_tr': {},
                '_cn': {},
                '_cr': {},
                '_tt': {},
                '_tc': {},
                '_ct': {},
                '_cc': {},
            }
        },
        'z': {
            'exts': {
                '_nn': {},
                '_nr': {},
                '_rn': {},
                '_rr': {},
                '_nt': {},
                '_nc': {},
                '_rt': {},
                '_rc': {},
                '_tn': {},
                '_tr': {},
                '_cn': {},
                '_cr': {},
                '_tt': {},
                '_tc': {},
                '_ct': {},
                '_cc': {},
            }
        },
    },
  },
  { 'base': '?gemm_small_kernel_b0',
    'modes': {
        's': {
            'exts': {
                '_nn': {},
                '_nt': {},
                '_tn': {},
                '_tt': {},
            }
        },
        'd': {
            'exts': {
                '_nn': {},
                '_nt': {},
                '_tn': {},
                '_tt': {},
            }
        },
        'c': {
            'exts': {
                '_nn': {},
                '_nr': {},
                '_rn': {},
                '_rr': {},
                '_nt': {},
                '_nc': {},
                '_rt': {},
                '_rc': {},
                '_tn': {},
                '_tr': {},
                '_cn': {},
                '_cr': {},
                '_tt': {},
                '_tc': {},
                '_ct': {},
                '_cc': {},
            }
        },
        'z': {
            'exts': {
                '_nn': {},
                '_nr': {},
                '_rn': {},
                '_rr': {},
                '_nt': {},
                '_nc': {},
                '_rt': {},
                '_rc': {},
                '_tn': {},
                '_tr': {},
                '_cn': {},
                '_cr': {},
                '_tt': {},
                '_tc': {},
                '_ct': {},
                '_cc': {},
            }
        },
    },
  },
]

search_order = [base_dict]

if conf_hdat.has('ARCH_X86_64')
  search_order = [x86_64_base_dict] + search_order
  if conf_hdat.has('HASWELL')
    search_order = [x86_64_haswell_dict] + search_order
  elif conf_hdat.has('SKYLAKEX')
    search_order = [x86_64_skylakex_dict, x86_64_haswell_dict] + search_order
  elif conf_hdat.has('ZEN')
    search_order = [x86_64_zen_dict] + search_order
  elif conf_hdat.has('SANDYBRIDGE')
    search_order = [x86_64_sandybridge_dict] + search_order
  endif
elif conf_hdat.has('ARCH_ARM64')
  search_order = [arm64_base_dict] + search_order
  if conf_hdat.has('ARMV8')
    search_order = [arm64_armv8_dict] + search_order
  endif
endif

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

      if base == '?gemm3m' and conf_hdat.has('ARCH_ARM64')
        continue
      endif

      src = ''
      foreach dict : search_order
        if dict.has_key(base) and dict[base].has_key(mode) and dict[base][mode].has_key(ext)
          src = dict[base][mode][ext]
          break
        endif
      endforeach
      if src == ''
        error(f'Missing src file for @base@ @mode@ @ext@', search_order, conf_hdat.keys())
      endif

      if extdat.has_key('addl')
        _ext_cargs += extdat['addl']
      endif
      sym_name = base.replace('?', mode) + ext
      sym_underscored = f'@sym_name@_'
      _ext_cargs += [
        f'-DASMNAME=@asm_name_prefix@@sym_name@',
        f'-DASMFNAME=@asm_name_prefix@@sym_underscored@',
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
