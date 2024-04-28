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
        'dir': 'x86_64',
        'kernel': 'swap.S',
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
        'dir': 'x86_64',
        'kernel': 'scal.S',
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
        'dir': 'x86_64',
        'kernel': 'axpy.S',
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
        'dir': 'x86_64',
        'kernel': 'dot.S',
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
  # # TODO: Handle the i*amax style
  # {'base': 'axpby',
  #  'modes': {'s': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  #            'd': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  #            'c': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  #            'z': {'dir': 'arm', 'kernel': 'axpby.c', 'exts': ['_k']},
  #           },
  # },
  # # Level 2 BLAS
  # # There are additional sources so now we have srcs
  # # Ordered as per KERNEL.generic and Makefile.L2
  # # exts are used to find the flags for each modality
  # # ext is attached to base (only useful from Level 2)
  # # i.e. baseext (e.g., gemv_n, gemv_t, cgeru_k, cgerc_k)
  # {'base': 'gemv',
  #  'modes': {'s': {'dir': 'arm', 'kernel': 'gemv_n.c', 'exts': ['_n', '_t']},
  #            'd': {'dir': 'arm', 'kernel': 'gemv_n.c',
  #                  # TODO: _t should take fmaflag, but then we need a dictionary..
  #                  # From Makefile.L2
  #                  'exts': ['_n', '_t']},
  #            'q': {'dir': 'arm', 'kernel': 'gemv_n.c', 'exts': ['_n', '_t']},
  #            'c': {'dir': 'arm', 'kernel': 'zgemv_n.c',
  #                  'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
  #            'z': {'dir': 'arm', 'kernel': 'zgemv_n.c',
  #                  'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
  #            'x': {'dir': 'arm', 'kernel': 'zgemv_n.c',
  #                  'exts': ['_n', '_t', '_r', '_c', '_o', '_u', '_s', '_d']},
  #           },
  # },
  # {'base': 'symv',
  #  'modes': {'s': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  #            'd': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  #            'q': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  #            'c': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  #            'z': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  #            'x': {'dir': 'generic', 'kernel': 'symv_k.c', 'exts': ['_U', '_L']},
  #           },
  # },
  # {'base': 'ger',
  #  'modes': {'s': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']},
  #            'd': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']},
  #            'q': {'dir': 'generic', 'kernel': 'ger.c', 'exts': ['_k']}
  #           },
  # },
  # {'base': 'geru',
  #  'modes': {'c': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  #            'z': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  #            'x': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']}
  #           }
  # },
  # {'base': 'gerc',
  #  'modes': {'c': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  #            'z': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  #            'x': {'dir': 'generic', 'kernel': 'zger.c', 'exts': ['_k']},
  #           },
  # },
  # {'base': 'hemv',
  #  'modes': {'c': {'dir': 'generic',
  #                  'kernel': 'zhemv_k.c',
  #                  'exts': ['_U', '_L', '_V', '_M']},
  #            'z': {'dir': 'generic',
  #                  'kernel': 'zhemv_k.c',
  #                  'exts': ['_U', '_L', '_V', '_M']},
  #            'x': {'dir': 'generic',
  #                  'kernel': 'zhemv_k.c',
  #                  'exts': ['_U', '_L', '_V', '_M']},
  #           },
  # },
  # {'base': 'bgemv',
  #  'modes': {'s': {'dir': 'x86_64',
  #                  'kernel': 'sbgemv_n.c',
  #                  'exts': ['_n', '_t']}
  #           },
  # },
]

kernel_confs = []
foreach _kop : base_kops
  base = _kop['base']
  modes = _kop['modes']
  foreach mode, details : modes
    __cargs = _cargs
    prec_mode = precision_mappings[mode]
    src = join_paths(details['dir'], details['kernel'])
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
    # Now the rest, one run for each ext
    foreach ext : details['exts']
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
