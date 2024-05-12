# Naming conventions can be read from here:
# https://click.rgoswami.me/intel_blas_names
# Ported in order from the Makefile
# TODO(rg): Add the CONJ flags
driver_kops = [
  { 'base': '?gbmv',
    'sources': {
        'gbmv_k.c': {'mode': ['s', 'd', 'q'], 'exts': ['_n', '_t']},
        'zgbmv_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_n', '_t', '_r', '_c',
                               '_o', '_u', '_s', '_d']},
    }
  },
  { 'base': '?gbmv_thread',
    'sources': {
        'gbmv_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_n', '_t', '_r', '_c',
                                   '_o', '_u', '_s', '_d']},
    }
  },
  { 'base': '?gemv_thread',
    'sources': {
        'gemv_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_n', '_t', '_r', '_c',
                                   '_o', '_u', '_s', '_d']},
    }
  },
  { 'base': '?ger_thread',
    'sources': {
        'ger_thread.c': {'mode': ['s', 'd', 'q'], 'exts': ['']},
    }
  },
  { 'base': '?ger_thread',
    'sources': {
        'ger_thread.c': {'mode': ['c', 'z', 'x'],
                         'exts': ['_U', '_C',
                                  '_V', '_D']},
    }
  },
  { 'base': '?symv_thread',
    'sources': {
        'symv_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_U', '_L']},
    }
  },
  { 'base': '?hemv_thread',
    'sources': {
        'symv_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_U', '_L',
                                   '_V', '_M']},
    },
  },
  { 'base': '?syr_thread',
    'sources': {
        'syr_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_U', '_L']},
    },
  },
  { 'base': '?her_thread',
    'sources': {
        'syr_thread.c': {'mode': ['c', 'z', 'x'],
                         'exts': ['_U', '_L',
                                  '_V', '_M']},
    },
  },
  { 'base': '?syr2_thread',
    'sources': {
        'syr2_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_U', '_L']},
    },
  },
  { 'base': '?her2_thread',
    'sources': {
        'syr2_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_U', '_L',
                                   '_V', '_M']},
    },
  },
  { 'base': '?hbmv',
    'sources': {
        'zhbmv_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L',
                               '_V', '_M']},
    },
  },
  { 'base': '?hbmv_thread',
    'sources': {
        'sbmv_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_U', '_L',
                                   '_V', '_M']},
    },
  },
  { 'base': '?her',
    'sources': {
        'zher_k.c': {'mode': ['c', 'z', 'x'],
                     'exts': ['_U', '_L',
                              '_V', '_M']},
    },
  },
  { 'base': '?her2',
    'sources': {
        'zher2_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L',
                               '_V', '_M']},
    },
  },
  { 'base': '?hpmv',
    'sources': {
        'zhpmv_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L',
                               '_V', '_M']},
    },
  },
  { 'base': '?hpmv_thread',
    'sources': {
        'spmv_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_U', '_L',
                                   '_V', '_M']},
    },
  },
  { 'base': '?hpr',
    'sources': {
        'zhpr_k.c': {'mode': ['c', 'z', 'x'],
                     'exts': ['_U', '_L',
                              '_V', '_M']},
    },
  },
  { 'base': '?hpr_thread',
    'sources': {
        'spr_thread.c': {'mode': ['c', 'z', 'x'],
                         'exts': ['_U', '_L',
                                  '_V', '_M']},
    },
  },
  { 'base': '?hpr2',
    'sources': {
        'zhpr2_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L',
                               '_V', '_M']},
    },
  },
  { 'base': '?hpr2_thread',
    'sources': {
        'spr2_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_U', '_L',
                                   '_V', '_M']},
    },
  },
  { 'base': '?sbmv',
    'sources': {
        'sbmv_k.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_U', '_L']},
        'zsbmv_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L']},
    }
  },
  { 'base': '?sbmv_thread',
    'sources': {
        'sbmv_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_U', '_L']},
    }
  },
  { 'base': '?spmv',
    'sources': {
        'spmv_k.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_U', '_L']},
        'zspmv_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L']},
    }
  },
  { 'base': '?spmv_thread',
    'sources': {
        'spmv_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_U', '_L']},
    }
  },
  { 'base': '?spr',
    'sources': {
        'spr_k.c': {'mode': ['s', 'd', 'q'],
                    'exts': ['_U', '_L']},
        'zspr_k.c': {'mode': ['c', 'z', 'x'],
                     'exts': ['_U', '_L']},
    }
  },
  { 'base': '?spr_thread',
    'sources': {
        'spr_thread.c': {'mode': ['s', 'd', 'q',
                                  'c', 'z', 'x'],
                         'exts': ['_U', '_L']},
    }
  },
  { 'base': '?spr2',
    'sources': {
        'spr2_k.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_U', '_L']},
        'zspr2_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L']},
    }
  },
  { 'base': '?spr2_thread',
    'sources': {
        'spr2_thread.c': {'mode': ['s', 'd', 'q',
                                   'c', 'z', 'x'],
                          'exts': ['_U', '_L']},
    }
  },
  { 'base': '?syr',
    'sources': {
        'syr_k.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_U', '_L']},
        'zsyr_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L']},
    }
  },
  { 'base': '?syr2',
    'sources': {
        'syr2_k.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_U', '_L']},
        'zsyr2_k.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_U', '_L']},
    }
  },
  { 'base': '?tbmv',
    'sources': {
        'tbmv_U.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NUU', '_NUN',
                              '_TLU', '_TLN']},
        'tbmv_L.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NLU', '_NLN',
                              '_TUU', '_TUN']},
        'ztbmv_U.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_NUU', '_NUN',
                               '_TLU', '_TLN',
                               '_CLU', '_CLN',
                               '_RUU', '_RUN']},
        'ztbmv_L.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_RLU', '_RLN',
                               '_NLU', '_NLN',
                               '_TUU', '_TUN',
                               '_CUU', '_CUN']},
    }
  },
  { 'base': '?tbmv_thread',
    'sources': {
        'tbmv_thread.c': {'mode': ['s', 'd', 'q'],
                          'exts': ['_NUU', '_NUN',
                                   '_TLU', '_TLN',
                                   '_NLU', '_NLN',
                                   '_TUU', '_TUN']},
    }
  },
  { 'base': '?tbmv_thread',
    'sources': {
        'tbmv_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_NUU', '_NUN',
                                   '_TLU', '_TLN',
                                   '_RLU', '_RLN',
                                   '_CLU', '_CLN',
                                   '_NLU', '_NLN',
                                   '_TUU', '_TUN',
                                   '_RUU', '_RUN',
                                   '_CUU', '_CUN']},
    }
  },
  { 'base': '?tbsv',
    'sources': {
        'tbsv_U.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NUU', '_NUN',
                              '_TLU', '_TLN']},
        'tbsv_L.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NLU', '_NLN',
                              '_TUU', '_TUN']},
        'ztbsv_U.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_NUU', '_NUN',
                               '_TLU', '_TLN',
                               '_CLU', '_CLN',
                               '_RUU', '_RUN']},
        'ztbsv_L.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_RLU', '_RLN',
                               '_NLU', '_NLN',
                               '_TUU', '_TUN',
                               '_CUU', '_CUN']},
    }
  },
  { 'base': '?tpmv',
    'sources': {
        'tpmv_U.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NUU', '_NUN',
                              '_TLU', '_TLN']},
        'tpmv_L.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NLU', '_NLN',
                              '_TUU', '_TUN']},
        'ztpmv_U.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_NUU', '_NUN',
                               '_TLU', '_TLN',
                               '_CLU', '_CLN',
                               '_RUU', '_RUN']},
        'ztpmv_L.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_RLU', '_RLN',
                               '_NLU', '_NLN',
                               '_TUU', '_TUN',
                               '_CUU', '_CUN']},
    }
  },
  { 'base': '?tpmv_thread',
    'sources': {
        'tpmv_thread.c': {'mode': ['s', 'd', 'q'],
                          'exts': ['_NUU', '_NUN',
                                   '_TLU', '_TLN',
                                   '_NLU', '_NLN',
                                   '_TUU', '_TUN']},
    }
  },
  { 'base': '?tpmv_thread',
    'sources': {
        'tpmv_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_NUU', '_NUN',
                                   '_TLU', '_TLN',
                                   '_RLU', '_RLN',
                                   '_CLU', '_CLN',
                                   '_NLU', '_NLN',
                                   '_TUU', '_TUN',
                                   '_RUU', '_RUN',
                                   '_CUU', '_CUN']},
    }
  },
  { 'base': '?tpsv',
    'sources': {
        'tpsv_U.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NUU', '_NUN',
                              '_TLU', '_TLN']},
        'tpsv_L.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NLU', '_NLN',
                              '_TUU', '_TUN']},
        'ztpsv_U.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_NUU', '_NUN',
                               '_TLU', '_TLN',
                               '_CLU', '_CLN',
                               '_RUU', '_RUN']},
        'ztpsv_L.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_RLU', '_RLN',
                               '_NLU', '_NLN',
                               '_TUU', '_TUN',
                               '_CUU', '_CUN']},
    }
  },
  { 'base': '?trmv',
    'sources': {
        'trmv_U.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NUU', '_NUN',
                              '_TLU', '_TLN']},
        'trmv_L.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NLU', '_NLN',
                              '_TUU', '_TUN']},
        'ztrmv_U.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_NUU', '_NUN',
                               '_TLU', '_TLN',
                               '_CLU', '_CLN',
                               '_RUU', '_RUN']},
        'ztrmv_L.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_RLU', '_RLN',
                               '_NLU', '_NLN',
                               '_TUU', '_TUN',
                               '_CUU', '_CUN']},
    }
  },
  { 'base': '?trmv_thread',
    'sources': {
        'trmv_thread.c': {'mode': ['s', 'd', 'q'],
                          'exts': ['_NUU', '_NUN',
                                   '_TLU', '_TLN',
                                   '_NLU', '_NLN',
                                   '_TUU', '_TUN']},
    }
  },
  { 'base': '?trmv_thread',
    'sources': {
        'trmv_thread.c': {'mode': ['c', 'z', 'x'],
                          'exts': ['_NUU', '_NUN',
                                   '_TLU', '_TLN',
                                   '_RLU', '_RLN',
                                   '_CLU', '_CLN',
                                   '_NLU', '_NLN',
                                   '_TUU', '_TUN',
                                   '_RUU', '_RUN',
                                   '_CUU', '_CUN']},
    }
  },
  { 'base': '?trsv',
    'sources': {
        'trsv_U.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NUU', '_NUN',
                              '_TLU', '_TLN']},
        'trsv_L.c': {'mode': ['s', 'd', 'q'],
                     'exts': ['_NLU', '_NLN',
                              '_TUU', '_TUN']},
        'ztrsv_U.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_NUU', '_NUN',
                               '_TLU', '_TLN',
                               '_CLU', '_CLN',
                               '_RUU', '_RUN']},
        'ztrsv_L.c': {'mode': ['c', 'z', 'x'],
                      'exts': ['_RLU', '_RLN',
                               '_NLU', '_NLN',
                               '_TUU', '_TUN',
                               '_CUU', '_CUN']},
    }
  },
  # TODO(rg): Add the bfloat conditionals from Makefile:3709
]

# Initialize kernel configurations list
kernel_confs = []

# Iterate through each kernel operation
foreach _kop : driver_kops
  base = _kop['base']
  sources = _kop['sources']

  # Generate the symbol flags
  _ckop_args = []
  if symb_defs.has_key(base)
    symb_base = symb_defs[base]
    if symb_base.has_key('def')
      foreach _d : symb_base['def']
        _ckop_args += ['-D' + _d]
      endforeach
    endif
    if symb_base.has_key('undef')
      foreach _u : symb_base['undef']
        _ckop_args += ['-U' + _u]
      endforeach
    endif
  endif

  # Iterate through each source file and its details
  foreach fname, details : sources
    modes = details['mode']
    exts = details['exts']

    # Iterate through each mode
    foreach mode : modes
      # Generate the mapping for the type
      __cargs = _cargs + _ckop_args
      prec_mode = precision_mappings[mode]

      # Add precision-specific definitions
      if prec_mode.has_key('def')
        foreach _d : prec_mode['def']
          __cargs += ['-D' + _d]
        endforeach
      endif
      if prec_mode.has_key('undef')
        foreach _u : prec_mode['undef']
          __cargs += ['-U' + _u]
        endforeach
      endif

      # Iterate through each extension
      foreach ext : exts
        _ext_cargs = []

        # Check ext_mappings first
        if ext_mappings.has_key(ext)
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
          foreach ext_map : ext_mappings_l2
            if ext_map['ext'] == ext and mode in ext_map['for']
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

        # Construct the final paths
        src = fname

        # Add additional flags if present
        if sources[fname].has_key('addl')
          _ext_cargs += sources[fname]['addl']
        endif

        # Generate the symbol name
        sym_name = base.replace('?', mode) + ext
        sym_underscored = f'@sym_name@_'

        # Add standard flags for naming conventions
        _ext_cargs += [
          f'-DASMNAME=@sym_name@',
          f'-DASMFNAME=@sym_underscored@',
          f'-DNAME=@sym_underscored@',
          f'-DCNAME=@sym_name@',
          f'-DCHAR_NAME="@sym_underscored@"',
          f'-DCHAR_CNAME="@sym_name@"',
        ]

        # Append the current configuration
        current_def = {
          'c_args': __cargs + _ext_cargs,
          'name': sym_name,
          'src': src,
        }
        kernel_confs += [current_def]
      endforeach
    endforeach
  endforeach
endforeach

# Create the static libraries from the configurations
_kern_libs = []
foreach conf : kernel_confs
    # message(conf['name'])
    # message(conf)
    _kern_libs += [static_library(
        conf['name'],
        conf['src'],
        include_directories: _inc,
        c_args: conf['c_args'],
    )]
endforeach

# Create the final kernel library
_l2_driver = static_library('l2_driver',
                           link_whole: _kern_libs)
