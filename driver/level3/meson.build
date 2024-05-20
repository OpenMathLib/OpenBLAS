# Naming conventions can be read from here:
# https://click.rgoswami.me/intel_blas_names
# Ported in order from the Makefile
driver_kops = [
  # { 'base': '?bgemm',
  #   'sources': {
  #       'gemm.c': {'mode': ['s'], 'srcs': ['level3.c'],
  #                  'exts': ['_nn', '_nt',
  #                           '_tn', '_tt']},
  #   }
  # },
  { 'base': '?gemm',
    'sources': {
        'gemm.c': {'mode': ['s', 'd',], #'q'],
                   'exts': ['_nn', '_nt',
                            '_tn', '_tt']},
    }
  },
  { 'base': '?gemm',
    'sources': {
        'gemm.c': {'mode': ['c', 'z',],# 'x'],
                   'srcs': ['level3.c'],
                   'exts': ['_nn', '_nt',
                            '_nr', '_nc',
                            '_tn', '_tt',
                            '_tr', '_tc',
                            '_rn', '_rt',
                            '_rr', '_rc',
                            '_cn', '_ct',
                            '_cr', '_cc']},
    }
  },
  { 'base': '?gemm_thread',
    'sources': {
        'gemm_thread_m.c': {'mode': [''],
                            'exts': ['_m']},
        'gemm_thread_n.c': {'mode': [''],
                            'exts': ['_n']},
        'gemm_thread_mn.c': {'mode': [''],
                            'exts': ['_mn']},
        'gemm_thread_variable.c': {'mode': [''],
                            'exts': ['_variable']},
    }
  },
  # { 'base': '?bgemm_thread',
  #   'sources': {
  #       'gemm.c': {'mode': ['s'], 'srcs': ['level3.c'],
  #                  'exts': ['_nn', '_nt',
  #                           '_tn', '_tt']},
  #   }
  # },
  { 'base': '?gemm_thread',
    'sources': {
        'gemm.c': {'mode': ['s', 'd'], #'q'],
                   'addl': ['-DTHREADED_LEVEL3'],
                   'exts': ['_nn', '_nt',
                            '_tn', '_tt']},
    }
  },
  { 'base': '?gemm_thread',
    'sources': {
        'gemm.c': {'mode': ['c', 'z'],# 'x'],
                   'addl': ['-DTHREADED_LEVEL3'],
                   'exts': ['_nn', '_nt',
                            '_nr', '_nc',
                            '_tn', '_tt',
                            '_tr', '_tc',
                            '_rn', '_rt',
                            '_rr', '_rc',
                            '_cn', '_ct',
                            '_cr', '_cc']},
    }
  },
  { 'base': '?trmm',
    'sources': {
        'trmm_L.c': {'mode': ['s', 'd'],# 'q'],
                     'exts': ['_LNUU', '_LNUN',
                              '_LNLU', '_LNLN',
                              '_LTUU', '_LTUN',
                              '_LTLU', '_LTLN']},
        'trmm_R.c': {'mode': ['s', 'd'],# 'q'],
                     'exts': ['_RNUU', '_RNUN',
                              '_RNLU', '_RNLN',
                              '_RTUU', '_RTUN',
                              '_RTLU', '_RTLN']},
    }
  },
  { 'base': '?trmm',
    'sources': {
        'trmm_L.c': {'mode': ['c', 'z'],# 'x'],
                     'exts': ['_LNUU', '_LNUN',
                              '_LNLU', '_LNLN',
                              '_LTUU', '_LTUN',
                              '_LTLU', '_LTLN',
                              '_LRUU', '_LRUN',
                              '_LRLU', '_LRLN',
                              '_LCUU', '_LCUN',
                              '_LCLU', '_LCLN',
                             ]},
        'trmm_R.c': {'mode': ['c', 'z'],# 'x'],
                     'exts': ['_RNUU', '_RNUN',
                              '_RNLU', '_RNLN',
                              '_RTUU', '_RTUN',
                              '_RTLU', '_RTLN',
                              '_RRUU', '_RRUN',
                              '_RRLU', '_RRLN',
                              '_RCUU', '_RCUN',
                              '_RCLU', '_RCLN',
                             ]},
    }
  },
  { 'base': '?symm',
    'sources': {
        'symm_k.c': {'mode': ['s', 'd', # 'q'],
                              'c', 'z'], # 'x'],
                     'exts': ['_LU', '_LL',
                              '_RU', '_RL']},
    }
  },
  { 'base': '?symm_thread',
    'sources': {
        'symm_k.c': {'mode': ['s', 'd', # 'q'],
                              'c', 'z'], # 'x'],
                     'addl': ['-DTHREADED_LEVEL3'],
                     'exts': ['_LU', '_LL',
                              '_RU', '_RL']},
    }
  },
  { 'base': '?syrk',
    'sources': {
        'syrk_k.c': {'mode': ['s', 'd', # 'q'],
                              'c', 'z'], # 'x'],
                     'exts': ['_UN', '_UT',
                              '_LN', '_LT']},
    }
  },
  { 'base': '?syrk_thread',
    'sources': {
        'syrk_k.c': {'mode': ['s', 'd', # 'q'],
                              'c', 'z'], # 'x'],
                     'addl': ['-DTHREADED_LEVEL3'],
                     'exts': ['_UN', '_UT',
                              '_LN', '_LT', '']},
    }
  },
  { 'base': '?syrk_kernel',
    'sources': {
        'syrk_kernel.c': {'mode': ['s', 'd', # 'q'],
                                   'c', 'z'], # 'x'],
                     'exts': ['_U', '_L']},
    }
  },
  { 'base': '?syr2k',
    'sources': {
        'syr2k_k.c': {'mode': ['s', 'd', # 'q'],
                               'c', 'z'], # 'x'],
                      'exts': ['_UN', '_UT',
                               '_LN', '_LT']},
    }
  },
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
        src = [ fname ]
        # if sources[fname].has_key('srcs')
        #   src += sources[fname]['srcs']
        # endif

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
    message(conf['name'])
    message(conf)
    _kern_libs += [static_library(
        conf['name'],
        conf['src'],
        include_directories: _inc,
        c_args: conf['c_args'],
    )]
endforeach

# Create the final kernel library
_l3_driver = static_library('l3_driver',
                           link_whole: _kern_libs)
