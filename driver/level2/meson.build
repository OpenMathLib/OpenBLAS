# Naming conventions can be read from here:
# https://click.rgoswami.me/intel_blas_names
# Ported in order from the Makefile
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
        'symv_thread.c': {'mode': ['c', 'z'],
                          'exts': ['_U', '_L',
                                   '_V', '_M']},
    },
  },
]
