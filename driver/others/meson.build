# Define the operations and their sources
others_ops = [
  { 'base': 'memory', 'source': 'memory.c' },
  { 'base': 'xerbla', 'source': 'xerbla.c' },
  { 'base': 'openblas_set_num_threads', 'source': 'openblas_set_num_threads.c' },
  { 'base': 'openblas_get_num_threads', 'source': 'openblas_get_num_threads.c', 'addl': [
    '-O2', '-DSMALL_MATRIX_OPT', '-DMAX_PARALLEL_NUMBER=1', '-DMAX_STACK_ALLOC=2048', '-Wall', '-DF_OTHERS_GFORT',
    '-fPIC', '-DSMP_SERVER', '-DNO_WARMUP', '-DMAX_CPU_NUMBER=12', '-DVERSION="0.3.26.dev"', '-msse3', '-mssse3',
    '-msse4.1', '-mavx', '-mavx2'
  ]},
  { 'base': 'openblas_get_num_procs', 'source': 'openblas_get_num_procs.c' },
  { 'base': 'openblas_get_parallel', 'source': 'openblas_get_parallel.c' },
  { 'base': 'openblas_error_handle', 'source': 'openblas_error_handle.c' },
  { 'base': 'openblas_env', 'source': 'openblas_env.c' },
  { 'base': 'blas_server', 'source': 'blas_server.c' },
  { 'base': 'divtable', 'source': 'divtable.c', 'addl': ['-UDOUBLE'] },
  { 'base': 'blasL1thread', 'source': 'blas_l1_thread.c' },
  { 'base': 'parameter', 'source': 'parameter.c' }
]

# Initialize configurations list
others_confs = []

# Iterate through each operation
foreach op : others_ops
  base = op['base']
  source = op['source']
  addl = op.get('addl', [])

  # Generate the symbol flags
  ckop_args = []
  if symb_defs.has_key(base)
    symb_base = symb_defs[base]
    if symb_base.has_key('def')
      foreach d : symb_base['def']
        ckop_args += ['-D' + d]
      endforeach
    endif
    if symb_base.has_key('undef')
      foreach u : symb_base['undef']
        ckop_args += ['-U' + u]
      endforeach
    endif
  endif

  # Default compilation arguments
  c_args = _cargs + ckop_args + addl + '-DMAX_PARALLEL_NUMBER=1'

  # Generate the symbol name
  sym_name = base
  sym_underscored = '@0@_'.format(sym_name)

  # Add standard flags for naming conventions
  c_args += [
    '-DASMNAME=@0@'.format(sym_name),
    '-DASMFNAME=@0@_'.format(sym_name),
    '-DNAME=@0@_'.format(sym_name),
    '-DCNAME=@0@'.format(sym_name),
    '-DCHAR_NAME="@0@_"'.format(sym_name),
    '-DCHAR_CNAME="@0@"'.format(sym_name)
  ]

  # Append the current configuration
  current_def = {
    'c_args': c_args,
    'name': sym_name,
    'src': [source]
  }
  others_confs += [current_def]
endforeach

# Create the static libraries from the configurations
others_libs = []
foreach conf : others_confs
  others_libs += [static_library(
    conf['name'],
    conf['src'],
    include_directories: _inc,
    c_args: conf['c_args']
  )]
endforeach

# Create the final interface library
_others = static_library('_interface', link_whole: others_libs)
