configurations = [
  # TODO: Also add defines for DOUBLE for d and the rest
  {'defs': [], 'name': 'imax'},  # Original, no defines
  {'defs': ['-DUSE_ABS'], 'name': 'iamax'},
  {'defs': ['-DUSE_ABS', '-DUSE_MIN'], 'name': 'iamin'},
  {'defs': ['-DUSE_MIN'], 'name': 'imin'},
]

_static_libs = []
foreach conf: configurations
    _static_libs += static_library(
        conf['name'],
        'imax.c',
        include_directories:  _inc,
        c_args: conf['defs'],
    )
endforeach

_interface = static_library('_interface',
                    link_whole: _static_libs)
