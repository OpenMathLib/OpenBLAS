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
  { 'base': 'rot',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'rot.c', 'addl': [fma3_flag]}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'rot.c', 'addl': [fma3_flag]}}},
        'q': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zrot.c'}}},
        'cs': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zrot.c'}}},
        'zd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zrot.c'}}},
        'xq': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zrot.c'}}},
    },
  },
  { 'base': 'swap',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'swap.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'swap.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'swap.S'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zswap.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zswap.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zswap.c'}}},
    },
  },
  { 'base': 'scal',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'scal.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'scal.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'scal.S'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zscal.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zscal.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zscal.c'}}},
    },
  },
  { 'base': 'copy',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'copy.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'copy.c'}}},
        'q': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zcopy.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zcopy.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zcopy.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zcopy.c'}}},
    },
  },
  { 'base': 'axpy',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'axpy.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'axpy.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'axpy.S'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zaxpy.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zaxpy.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zaxpy.c'}}},
    },
  },
  { 'base': 'dot',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'dot.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'dot.c'}}},
        'q': {'exts': {'_k': {'dir': 'x86_64', 'kernel': 'dot.S'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zdot.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zdot.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zdot.c'}}},
    },
  },
  # xDOTU xDOTC xxDOT aren't present
  { 'base': 'nrm2',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'nrm2.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'nrm2.c'}}},
        'q': {'exts': {'_k': {'dir': 'arm', 'kernel': 'znrm2.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'znrm2.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'znrm2.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'znrm2.c'}}},
    },
  },
  { 'base': 'asum',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'asum.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'asum.c'}}},
        'q': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zasum.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zasum.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zasum.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zasum.c'}}},
    },
  },
  { 'base': 'amax',
    'modes': {
        's': {'exts': {'_k': {'dir': 'arm', 'kernel': 'amax.c'}}},
        'd': {'exts': {'_k': {'dir': 'arm', 'kernel': 'amax.c'}}},
        'q': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamax.c'}}},
        'c': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamax.c'}}},
        'z': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamax.c'}}},
        'x': {'exts': {'_k': {'dir': 'arm', 'kernel': 'zamax.c'}}},
    },
  },
  # # # TODO: Handle the i*amax style
  # # {'base': 'axpby',
  # #  'modes': {'s': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  # #            'd': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  # #            'c': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  # #            'z': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  # #           },
  # # },
  # # # Level 2 BLAS
  # # # There are additional sources so now we have srcs
  # # # Ordered as per KERNEL.generic and Makefile.L2
  # # # exts are used to find the flags for each modality
  # # # ext is attached to base (only useful from Level 2)
  # # # i.e. baseext (e.g., gemv_n, gemv_t, cgeru_k, cgerc_k)
  # # {'base': 'gemv',
  # #  'modes': {'s': {'dir': 'arm', 'kernel': 'gemv_n.c', 'exts': ['_n', '_t']},
  # #            'd': {'dir': 'arm', 'kernel': 'gemv_n.c',
  # #                  # TODO: _t should take fmaflag, but then we need a dictionary..
  # #                  # From Makefile.L2
  # #                  'exts': ['_n', '_t']},
  # #            'q': {'dir': 'arm', 'kernel': 'gemv_n.c', 'exts': ['_n', '_t']},
  # #            'c': {'dir': 'arm', 'kernel': 'zgemv_n.c',
  # #                  'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
  # #            'z': {'dir': 'arm', 'kernel': 'zgemv_n.c',
  # #                  'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
  # #            'x': {'dir': 'arm', 'kernel': 'zgemv_n.c',
  # #                  'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
  # #           },
  # # },
  # # {'base': 'symv',
  # #  'modes': {'s': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  # #            'd': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  # #            'q': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  # #            'c': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  # #            'z': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  # #            'x': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  # #           },
  # # },
  # # {'base': 'ger',
  # #  'modes': {'s': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']},
  # #            'd': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']},
  # #            'q': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']}
  # #           },
  # # },
  # # {'base': 'geru',
  # #  'modes': {'c': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  # #            'z': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  # #            'x': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']}
  # #           }
  # # },
  # # {'base': 'gerc',
  # #  'modes': {'c': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  # #            'z': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  # #            'x': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  # #           },
  # # },
  # # {'base': 'hemv',
  # #  'modes': {'c': {'dir': 'generic',
  # #                  'kernel': 'zhemv_k.c',
  # #                  'exts': ['_U', '_L', '_V', '_M']},
  # #            'z': {'dir': 'generic',
  # #                  'kernel': 'zhemv_k.c',
  # #                  'exts': ['_U', '_L', '_V', '_M']},
  # #            'x': {'dir': 'generic',
  # #                  'kernel': 'zhemv_k.c',
  # #                  'exts': ['_U', '_L', '_V', '_M']},
  # #           },
  # # },
  # # {'base': 'bgemv',
  # #  'modes': {'s': {'dir': 'x86_64',
  # #                  'kernel': 'sbgemv_n.c',
  # #                  'exts': ['_n', '_t']}
  # #           },
  # # },
]

kernel_confs = []
foreach _kop : base_kops
  base = _kop['base']
  modes = _kop['modes']
  foreach mode, details : modes
    __cargs = _cargs
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
      extmap = ext_mappings[ext]
      if extmap.has_key('def')
        foreach _d : extmap['def']
          __cargs += ('-D' + _d)
        endforeach
      endif
      if extmap.has_key('undef')
        foreach _u : extmap['undef']
          __cargs += ('-U' + _u)
        endforeach
      endif
      # Construct the final paths
      src = join_paths(extdat['dir'], extdat['kernel'])
      if extdat.has_key('addl')
        __cargs += extdat['addl']
      endif
      message(mode + base + ext)
      sym_name = mode + base + ext
      sym_underscored = f'@sym_name@_'
      __cargs += [
        f'-DASMNAME=@sym_name@',
        f'-DASMFNAME=@sym_underscored@',
        f'-DNAME=@sym_underscored@',
        f'-DCNAME=@sym_name@',
        f'-DCHAR_NAME="@sym_underscored@"',
        f'-DCHAR_CNAME="@sym_name@"',
      ]
      current_def = {
        'c_args': __cargs,
        'name': sym_name,
        'src': src
      }
      message(__cargs)
      kernel_confs += current_def
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
