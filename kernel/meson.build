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
  # TODO(rg): Handle the if defines set in arch
  # { 'base': '?gemm_beta', # This is a bfloat16 symbol, skipping for now
  #   'modes': {
  #       's': {'exts': {'': {'dir': 'x86_64', 'kernel': 'gemm_beta.S'}}},
  #       'd': {'exts': {'': {'dir': 'x86_64', 'kernel': 'gemm_beta.S'}}},
  #       'c': {'exts': {'': {'dir': 'x86_64', 'kernel': 'zgemm_beta.S'}}},
  #       'z': {'exts': {'': {'dir': 'x86_64', 'kernel': 'zgemm_beta.S'}}},
  #       # 'q': {'exts': {'': {'dir': 'generic', 'kernel': 'gemm_beta.c'}}},
  #       # 'x': {'exts': {'': {'dir': 'generic', 'kernel': 'zgemm_beta.c'}}},
  #   },
  # },
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
      extmap = ext_mappings[ext]
      if extmap.has_key('def')
        foreach _d : extmap['def']
          _ext_cargs += ('-D' + _d)
        endforeach
      endif
      if extmap.has_key('undef')
        foreach _u : extmap['undef']
          _ext_cargs += ('-U' + _u)
        endforeach
      endif
      # Construct the final paths
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
