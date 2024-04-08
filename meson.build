# Conventions:
# _ implies that the variables are not meant to be used outside here
# Optionals are applied from the top-level meson_options.txt
# They are declared at the top
# Typically derived from (in order) the CMakeLists.txt and Makefiles
#
# Installation:
# meson setup build --buildtype release
# meson compile -C build
# meson install --prefix=$HOME/.local/lapack
#
# NOTE: This is still a work in progress, the Makefiles are canonical
project('OpenBLAS', 'c',
        default_options: ['c_std=c99'],
        version: '0.3.26.dev')

openblas_major_version = 0 # soversion
openblas_minor_version = 3
openblas_patch_version = '26.dev'
openblas_version = f'@openblas_major_version@.@openblas_minor_version@.@openblas_patch_version@'

# Skip the check for valid CC
cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

# Common args
_args = []
_cargs = []
_fargs = []

# System configuration
build_single = get_option('build_single')
build_double = get_option('build_double')
build_complex = get_option('build_complex')
build_complex16 = get_option('build_complex16')

# Options from CMakelists
build_without_lapack = get_option('build_without_lapack')
build_lapack_deprecated = get_option('build_lapack_deprecated')
build_testing = get_option('build_testing')
use_c_lapack = get_option('use_c_lapack')
build_without_cblas = get_option('build_without_cblas')
dynamic_arch = get_option('dynamic_arch')
dynamic_older = get_option('dynamic_older')
build_relapack = get_option('build_relapack')
use_locking = get_option('use_locking')
use_perl = get_option('use_perl')
no_warmup = get_option('no_warmup')
no_affinity = get_option('no_affinity')
build_cpp_thread_safety_test = get_option('build_cpp_thread_safety_test')
build_cpp_thread_safety_gemv = get_option('build_cpp_thread_safety_gemv')
build_static_libs = get_option('build_static_libs')
max_stack_alloc = get_option('max_stack_alloc')

if host_machine.system() == 'linux'
    no_affinity = true
else
    no_affinity = false
endif

# Makefile.prebuild stuff
# TODO: Handle cpuidemu, and the target flags
# getarch = executable('getarch',
#           ['getarch.c', 'cpuid.S'])
# getarch_two = executable('getarch_2nd',
#                          ['getarch_2nd.c'])
# config_h = custom_target('gen_config_h',
#                           # input: ['getarch.c'],
#                           output: ['config.h'],
#                           command: [getarch, '3']
#                           )

_check_prefix = []
conf_data = configuration_data()
is_win = build_machine.system() == 'win'
conf_data.set('OS_WINDOWS', is_win)

getarch_array = [
 # Redundant, we have build_machine.system() for windows
 {'checks': ['WIN32', 'WIN64', 'CYGWIN32'], 'addl': ['_WIN64', '_WIN32'], 'name': 'Windows', 'def': ['OS_WINDOWS']},
 {'checks': ['i386', 'x86_64'], 'addl': ['_M_IX86', '_M_X64'], 'name': 'Intel AMD', 'def': ['INTEL_AMD']},
 {'checks': ['powerpc', 'ppc', 'POWERPC'], 'addl': ['__powerpc', 'powerpc', 'ppc', '_POWER'], 'name': 'PowerPC', 'def': ['POWER', 'ARCH_POWER']},
 {'checks': ['zarch', 's390x'], 'name': 'zarch', 'def': ['ARCH_ZARCH', 'ZARCH']},
 {'checks': ['ia64'], 'name': 'IA64', 'def': ['ARCH_IA64']},
 {'name': 'Alpha', 'addl': ['__alpha'], 'def': ['ARCH_ALPHA']},
 {'name': 'SPARC', 'addl': ['sparc'], 'def': ['ARCH_SPARC']},
 {'checks': ['mips'], 'name': 'MIPS', 'def': ['ARCH_MIPS']},
 {'name': 'MIPS 64', 'addl': ['__mips64'], 'def': ['ARCH_MIPS64']},
 {'name': 'LOONGARCH 64', 'addl': ['__loongarch64'], 'def': ['ARCH_LOONGARCH64']},
 {'name': 'RISCV 64', 'addl': ['__riscv'], 'def': ['ARCH_RISCV64']},
 {'checks': ['arm'], 'name': 'ARM', 'def': ['ARCH_ARM']},
 {'checks': ['aarch64'], 'name': 'ARM 64', 'def': ['ARCH_ARM64']},
]

foreach arch : getarch_array
  check_str = ''
  if arch.has_key('checks')
    foreach check : arch['checks']
      if check_str == ''
        check_str += 'defined(__' + check + '__)'
      else
        check_str += ' || defined(__' + check + '__)'
      endif
    endforeach
  endif
  if arch.has_key('addl')
    foreach check : arch['addl']
      if check_str == ''
        check_str += 'defined(' + check + ')'
      else
        check_str += ' || defined(' + check + ')'
      endif
    endforeach
  endif

  sname = arch['name']
  code = f'''
    #if !(@check_str@)
    #error "Not @sname@"
    #endif
  '''
  # message(code)

  is_arch = cc.compiles(code, name: arch['name'])
  foreach def : arch['def']
    conf_data.set(def, is_arch)
  endforeach
endforeach

configure_file(output : 'getarch_conf.h',
               configuration : conf_data)

# run_target('generate_config_h',
#            command: [meson.current_build_dir() + '/getarch', '1'],
#            depends: getarch)

# gch = run_command(meson.current_build_dir() + '/getarch', '1')
# outp = gch.stdout().strip()

# conf_data = configuration_data()
# configure_file(input : meson.current_build_dir() + '/config.h',
#   output : 'config.h',
#   configuration : conf_data)

# Makefile.system
cpu_fam = target_machine.cpu_family()

if cpu_fam in ['x86_64', 'ppc64', 'ppc']
  add_project_arguments([
   '-DSMALL_MATRIX_OPT',
   f'-DMAX_STACK_ALLOC=@max_stack_alloc@',
  ], language:'c')
endif

if fc_id == 'gcc'
    add_project_arguments('-DF_INTERFACE_GFORT', language: 'c')
elif fc_id == 'g95'
    add_project_arguments('-DF_INTERFACE_G95', language: 'c')
elif fc_id == 'flang'
    add_project_arguments('-DF_INTERFACE_FLANG', language: 'c')
elif fc_id in ['intel', 'intel-cl']
    add_project_arguments('-DF_INTERFACE_INTEL', language: 'c')
elif fc_id == 'pathscale'
    add_project_arguments('-DF_INTERFACE_PATHSCALE', language: 'c')
elif fc_id == 'pgi'
    add_project_arguments('-DF_INTERFACE_PGI', language: 'c')
elif fc_id == 'sun'
    add_project_arguments('-DF_INTERFACE_SUN', language: 'c')
endif

# Ignoring other hostarch checks and conflicts for arch in BSD for now
_inc = include_directories('.')
# subdir('lapack-netlib')
subdir('interface')
subdir('kernel')

_openblas = static_library('openblas',
                    link_whole: [ _interface, _kern])
