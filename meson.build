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
project('OpenBLAS', ['c', 'fortran'],
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

# TODO: Be conditional
_cargs = [
 '-DBUILD_SINGLE=1',
 '-DBUILD_DOUBLE=1',
 '-DBUILD_COMPLEX=1',
 '-DBUILD_COMPLEX16=1',
]

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
is_win = host_machine.system() == 'windows' or host_machine.system() == 'cygwin'
conf_data.set('OS_WINDOWS', is_win)
hostcpu = host_machine.cpu_family()
conf_data.set('ARCH', hostcpu.to_upper())

defarch_array = [
 # {'system': ['windows', 'cygwin'], 'def': ['OS_WINDOWS']},
 {'cpu': ['aarch64'], 'def': ['ARCH_ARGM64']},
 {'cpu': ['alpha'], 'def': ['ARCH_ALPHA']},
 {'cpu': ['arm'], 'def': ['ARCH_ARM', 'USE_TRMM']},
 {'cpu': ['x86_64'], 'def': ['INTEL_AMD', 'USE_GEMM3M', 'USE_DIRECT_SGEMM']},
 {'cpu': ['s390x'], 'def': ['ARCH_ZARCH', 'ZARCH']},
 {'cpu': ['ia64'], 'def': ['ARCH_IA64', 'USE_GEMM3M']},
 {'cpu': ['sparc'], 'def': ['ARCH_SPARC']},
 {'cpu': ['mips'], 'def': ['ARCH_MIPS']},
 {'cpu': ['mips64'], 'def': ['ARCH_MIPS64', 'USE_TRMM']},
 {'cpu': ['loongarch64'], 'def': ['ARCH_LOONGARCH64']},
 {'cpu': ['riscv64'], 'def': ['ARCH_RISCV64', 'USE_TRMM']},
 {'cpu': ['ppc64', 'ppc'], 'def': ['ARCH_POWER', 'POWER']},
 {'cpu': ['generic'], 'def': ['USE_TRMM']},
]

# TODO(rg): Handle the kenel architectures later

foreach arch : defarch_array
  is_cpu = hostcpu in arch['cpu']
  foreach def : arch['def']
    conf_data.set(def, is_cpu)
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

if cpu_fam == 'x86_64'
  _cargs += ['-m64']
endif

if fc_id == 'gcc'
    add_project_arguments('-DF_INTERFACE_GFORT', language: 'c')
else
  upper_fcid = fc_id.to_upper()
  add_project_arguments(f'-DF_INTERFACE_@upper_fcid@', language: 'c')
endif

simd_extensions = [
  {'flag': '-mmmx', 'define': 'HAVE_MMX'},
  {'flag': '-msse', 'define': 'HAVE_SSE'},
  {'flag': '-msse2', 'define': 'HAVE_SSE2'},
  {'flag': '-msse3', 'define': 'HAVE_SSE3'},
  {'flag': '-mssse3', 'define': 'HAVE_SSSE3'},
  {'flag': '-msse4.1', 'define': 'HAVE_SSE4_1'},
  {'flag': '-msse4.2', 'define': 'HAVE_SSE4_2'},
  {'flag': '-mfma', 'define': 'HAVE_FMA3'},
  {'flag': '-mavx', 'define': 'HAVE_AVX'},
  {'flag': '-mavx2', 'define': 'HAVE_AVX2'},
]

simd_cargs = []
simd_conf_data = configuration_data()

foreach ext : simd_extensions
  if cc.has_argument(ext['flag'])
    simd_cargs += ext['flag']
    simd_conf_data.set(ext['define'], true)
  else
    simd_conf_data.set(ext['define'], false)
  endif
endforeach

# Generate configuration header
configure_file(output : 'simd_conf.h',
               configuration : simd_conf_data)
# Can't be added as a project argument since it won't be generated until later
bldroot = meson.build_root()
_cargs += ['-include', f'@bldroot@/simd_conf.h']
add_project_arguments(simd_cargs, language: 'c')

# Common symbol related options
symnames = ['ASMNAME', 'ASMFNAME', 'NAME', 'CNAME', 'CHAR_NAME', 'CHAR_CNAME']

# TODO(rg): Maybe make these conditional..
_cargs += [
 '-DSMP_SERVER', # This is evidently necessary for the driver/level2
 '-DBUILD_SINGLE=1',
 '-DBUILD_DOUBLE=1',
 '-DBUILD_COMPLEX=1',
 '-DBUILD_COMPLEX16=1',
]

# Other common options, move later
# Undefine to help prevent clashes
foreach symb : symnames
  _cargs += f'-U@symb@'
endforeach
# Based on options
if not build_without_cblas
  _cargs += '-DCBLAS'
endif
if no_affinity
  _cargs += '-DNO_AFFINITY'
endif
if no_warmup
  _cargs += '-DNO_WARMUP'
endif
# Parallel builds
# TODO: This can be cleaned up significantly
# Also use multiprocessing.cpu_count()
# TODO: Handle SMP_SERVER
num_parallel = get_option('num_parallel')
if num_parallel > 1
  _cargs += f'-DMAX_PARALLEL_NUMBER=@num_parallel@'
endif
num_cores = get_option('num_cores')
if num_cores > 0
  num_threads = num_cores
else
  num_threads = 0
endif
use_thread = false
if num_threads > 2
  use_thread = true
endif
if use_thread
  message('Multi-threading enabled with ' + num_threads.to_string() + ' threads.')
  _cargs += f'-DMAX_CPU_NUMBER=@num_threads@'
  _cargs += f'-DSMP_SERVER'
else
  if get_option('use_locking')
    _cargs += '-DUSE_LOCKING'
  endif
endif

# Common maps
# Naming conventions: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2024-1/naming-conventions-for-blas-routines.html
# Also see:
# L1: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2024-1/blas-level-1-routines-and-functions.html
precision_mappings = {
    's': {'undef': ['COMPLEX', 'DOUBLE']},
    'd': {'undef': ['COMPLEX'], 'def': ['DOUBLE']},
    'q': {'undef': ['COMPLEX'], 'def': ['XDOUBLE']},
    'c': {'undef': ['DOUBLE'], 'def': ['COMPLEX']},
    'z': {'def': ['COMPLEX', 'DOUBLE']},
    'x': {'def': ['COMPLEX', 'XDOUBLE']},
    'cs': {'undef': ['DOUBLE'], 'def': ['COMPLEX']},
    'sc': {'undef': ['DOUBLE'], 'def': ['COMPLEX']},
    'dz': {'def': ['COMPLEX', 'DOUBLE']},
    'zd': {'def': ['COMPLEX', 'DOUBLE']},
    'qx': {'def': ['COMPLEX', 'XDOUBLE']},
    'xq': {'def': ['COMPLEX', 'XDOUBLE']},
    '': {}, # special case, for cblas_?dot*_sub
    # NOTE: Anything with XDOUBLE aka longdouble has no cblas_
    # xq / qx == x
    # sc / cs == c
    # zd / dz is the same as z
    # 'zd': {'undef': [], 'def': ['COMPLEX', 'DOUBLE']},
}

ext_mappings = {
    # TRANSA is only in drivers level2, kernel uses TRANS
    # LOWER is not there for only for ?ger_thread_U in drivers level2
    # TODO(rg): Does that, i.e. having (un)used symbols (un)defined matter?
    '': {}, # special case
    '_k': {},
    '_U': {'undef': ['LOWER', 'CONJ', 'XCONJ']},
    '_C': {'def': ['CONJ'], 'undef': ['XCONJ']},
    # '_V': {'def': ['XCONJ'], 'undef': ['CONJ']},
    '_D': {'def': ['CONJ', 'XCONJ']},
    '_L': {'def': ['LOWER']},
    '_LN': {'def': ['LEFT'], 'undef': ['TRANSA'], 'except': ['?syrk']},
    # Handle HEMV and HEMVREV better
    '_V': {'def': ['HEMV', 'HEMVREV', 'XCONJ'], 'undef': ['LOWER', 'CONJ']},
    '_M': {'def': ['HEMV', 'HEMVREV', 'LOWER']},
    '_n': {'undef': ['TRANS', 'TRANSA', 'CONJ', 'XCONJ']},
    '_t': {'def': ['TRANS', 'TRANSA'], 'undef': ['CONJ', 'XCONJ']},
    '_r': {'def': ['CONJ'], 'undef': ['TRANS', 'TRANSA', 'XCONJ']},
    '_c': {'def': ['TRANS', 'TRANSA', 'CONJ'], 'undef': ['XCONJ']},
    '_o': {'def': ['XCONJ'], 'undef': ['TRANS', 'TRANSA', 'CONJ']},
    '_u': {'def': ['TRANS', 'TRANSA', 'XCONJ'], 'undef': ['CONJ']},
    '_s': {'def': ['CONJ', 'XCONJ'], 'undef': ['TRANS', 'TRANSA']},
    '_d': {'def': ['TRANS', 'TRANSA', 'CONJ', 'XCONJ']},
    '_nn': {'def': ['NN']},
    '_nt': {'def': ['NT']},
    '_nr': {'def': ['NR']},
    '_nc': {'def': ['NC']},
    '_tn': {'def': ['TN']},
    '_tt': {'def': ['TT']},
    '_tr': {'def': ['TR']},
    '_tc': {'def': ['TC']},
    '_rn': {'def': ['RN']},
    '_rt': {'def': ['RT']},
    '_rr': {'def': ['RR']},
    '_rc': {'def': ['RC']},
    '_cn': {'def': ['CN']},
    '_ct': {'def': ['CT']},
    '_cr': {'def': ['CR']},
    '_cc': {'def': ['CC']},
    # Level 3 symbols
    '_LU': {'def': ['NN'], 'undef': ['LOWER', 'RSIDE']},
    '_LL': {'def': ['LOWER', 'NN'], 'undef': ['RSIDE']},
    '_RU': {'def': ['RSIDE', 'NN'], 'undef': ['LOWER']},
    '_RL': {'def': ['RSIDE', 'NN', 'LOWER']},
    '_UN': {'undef': ['TRANS', 'LOWER'], 'except': ['?syrk']},
    '_UT': {'def': ['TRANS'], 'undef': ['LOWER'], 'except': ['?syrk']},
}

ext_mappings_l2 = [
    {'ext': '_NUU', 'def': ['UNIT'], 'undef': ['TRANSA'], 'for': ['s', 'd']},
    {'ext': '_NUN', 'undef': ['TRANSA', 'UNIT'], 'for': ['s', 'd']},
    {'ext': '_TLU', 'def': ['UNIT', 'TRANSA'], 'for': ['s', 'd']},
    {'ext': '_TLN', 'def': ['TRANSA'], 'undef': ['UNIT'], 'for': ['s', 'd']},
    {'ext': '_NLU', 'def': ['UNIT'], 'undef': ['TRANSA'], 'for': ['s', 'd']},
    {'ext': '_NLN', 'undef': ['TRANSA', 'UNIT'], 'for': ['s', 'd']},
    {'ext': '_TUU', 'def': ['UNIT', 'TRANSA'], 'for': ['s', 'd']},
    {'ext': '_TUN', 'def': ['TRANSA'], 'undef': ['UNIT'], 'for': ['s', 'd']},
    {'ext': '_NUU', 'def': ['UNIT', 'TRANSA=1'], 'for': ['c', 'x', 'z']},
    {'ext': '_NUN', 'def': ['TRANSA=1'], 'undef': ['UNIT'], 'for': ['c', 'x', 'z']},
    {'ext': '_TLU', 'def': ['UNIT', 'TRANSA=2'], 'for': ['c', 'x', 'z']},
    {'ext': '_TLN', 'def': ['TRANSA=2'], 'undef': ['UNIT'], 'for': ['c', 'x', 'z']},
    {'ext': '_RLU', 'def': ['UNIT', 'TRANSA=3'], 'for': ['c', 'x', 'z']},
    {'ext': '_RLN', 'def': ['TRANSA=3'], 'undef': ['UNIT'], 'for': ['c', 'x', 'z']},
    {'ext': '_CLU', 'def': ['UNIT', 'TRANSA=4'], 'for': ['c', 'x', 'z']},
    {'ext': '_CLN', 'def': ['TRANSA=4'], 'undef': ['UNIT'], 'for': ['c', 'x', 'z']},
    {'ext': '_NLU', 'def': ['UNIT', 'TRANSA=1'], 'for': ['c', 'x', 'z']},
    {'ext': '_NLN', 'def': ['TRANSA=1'], 'undef': ['UNIT'], 'for': ['c', 'x', 'z']},
    {'ext': '_TUU', 'def': ['UNIT', 'TRANSA=2'], 'for': ['s', 'd']},
    {'ext': '_TUN', 'def': ['TRANSA=2'], 'undef': ['UNIT'], 'for': ['s', 'd']},
    {'ext': '_RUU', 'def': ['UNIT', 'TRANSA=3'], 'for': ['c', 'x', 'z']},
    {'ext': '_RUN', 'def': ['TRANSA=3'], 'undef': ['UNIT'], 'for': ['c', 'x', 'z']},
    {'ext': '_CUU', 'def': ['UNIT', 'TRANSA=4'], 'for': ['c', 'x', 'z']},
    {'ext': '_CUN', 'def': ['TRANSA=4'], 'undef': ['UNIT'], 'for': ['c', 'x', 'z']}
]

ext_mappings_l3 = [
 {'ext': '_LNUU', 'def': ['UPPER', 'UNIT'], 'undef': ['TRANSA'], 'for': ['s', 'd']},
 {'ext': '_LNUN', 'def': ['UPPER'], 'undef': ['TRANSA', 'UNIT'], 'for': ['s', 'd']},
 {'ext': '_LNLU', 'def': ['UNIT'], 'undef': ['TRANSA', 'UPPER'], 'for': ['s', 'd']},
 {'ext': '_LNLN', 'undef': ['TRANSA', 'UPPER', 'UNIT'], 'for': ['s', 'd', 'q']},
 {'ext': '_LTUU', 'def': ['TRANSA', 'UPPER', 'UNIT'], 'for': ['s', 'd']},
 {'ext': '_LTUN', 'def': ['TRANSA', 'UPPER'], 'undef': ['UNIT'], 'for': ['s', 'd']},
 {'ext': '_LTLU', 'def': ['TRANSA', 'UNIT'], 'undef': ['UPPER'], 'for': ['s', 'd', 'q']},
 {'ext': '_LTLN', 'def': ['TRANSA'], 'undef': ['UPPER', 'UNIT'], 'for': ['s', 'd']},
 {'ext': '_RNUU', 'def': ['UPPER', 'UNIT'], 'undef': ['TRANSA'], 'for': ['s', 'd']},
 {'ext': '_RNUN', 'def': ['UPPER'], 'undef': ['TRANSA', 'UNIT'], 'for': ['s', 'd']},
 {'ext': '_RNLU', 'def': ['UNIT'], 'undef': ['TRANSA', 'UPPER'], 'for': ['s', 'd']},
 {'ext': '_RNLN', 'undef': ['TRANSA', 'UPPER', 'UNIT'], 'for': ['s', 'd']},
 {'ext': '_RTUU', 'def': ['TRANSA', 'UPPER', 'UNIT'], 'for': ['s', 'd']},
 {'ext': '_RTUN', 'def': ['TRANSA', 'UPPER'], 'undef': ['UNIT'], 'for': ['s', 'd', 'q']},
 {'ext': '_RTLU', 'def': ['TRANSA', 'UNIT'], 'undef': ['UPPER'], 'for': ['s', 'd', 'q']},
 {'ext': '_RTLN', 'def': ['TRANSA'], 'undef': ['UPPER', 'UNIT'], 'for': ['s', 'd']},

 # For trmm
 {'ext': '_LNUU', 'def': ['UPPER', 'UNIT'], 'undef': ['TRANSA', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LNUN', 'def': ['UPPER'], 'undef': ['TRANSA', 'UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LNLU', 'def': ['UNIT'], 'undef': ['TRANSA', 'UPPER', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LNLN', 'def': [], 'undef': ['TRANSA', 'UPPER', 'UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LTUU', 'def': ['TRANSA', 'UPPER', 'UNIT'], 'undef': ['CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LTUN', 'def': ['TRANSA', 'UPPER'], 'undef': ['UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LTLU', 'def': ['TRANSA', 'UNIT'], 'undef': ['UPPER', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LTLN', 'def': ['TRANSA'], 'undef': ['UPPER', 'UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_LRUU', 'def': ['UPPER', 'UNIT', 'CONJ'], 'undef': ['TRANSA'], 'for': ['c', 'z', 'x']},
 {'ext': '_LRUN', 'def': ['UPPER', 'CONJ'], 'undef': ['TRANSA', 'UNIT'], 'for': ['c', 'z', 'x']},
 {'ext': '_LRLU', 'def': ['UNIT', 'CONJ'], 'undef': ['TRANSA', 'UPPER'], 'for': ['c', 'z', 'x']},
 {'ext': '_LRLN', 'def': ['CONJ'], 'undef': ['TRANSA', 'UPPER', 'UNIT'], 'for': ['c', 'z', 'x']},
 {'ext': '_LCUU', 'def': ['TRANSA', 'UPPER', 'UNIT', 'CONJ'], 'undef': [], 'for': ['c', 'z', 'x']},
 {'ext': '_LCUN', 'def': ['TRANSA', 'UPPER', 'CONJ'], 'undef': ['UNIT'], 'for': ['c', 'z', 'x']},
 {'ext': '_LCLU', 'def': ['TRANSA', 'UNIT', 'CONJ'], 'undef': ['UPPER'], 'for': ['c', 'z', 'x']},
 {'ext': '_LCLN', 'def': ['TRANSA', 'CONJ'], 'undef': ['UPPER', 'UNIT'], 'for': ['c', 'z', 'x']},
 {'ext': '_RNUU', 'def': ['UPPER', 'UNIT'], 'undef': ['TRANSA', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RNUN', 'def': ['UPPER'], 'undef': ['TRANSA', 'UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RNLU', 'def': ['UNIT'], 'undef': ['TRANSA', 'UPPER', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RNLN', 'def': [], 'undef': ['TRANSA', 'UPPER', 'UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RTUU', 'def': ['TRANSA', 'UPPER', 'UNIT'], 'undef': ['CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RTUN', 'def': ['TRANSA', 'UPPER'], 'undef': ['UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RTLU', 'def': ['TRANSA', 'UNIT'], 'undef': ['UPPER', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RTLN', 'def': ['TRANSA'], 'undef': ['UPPER', 'UNIT', 'CONJ'], 'for': ['c', 'z', 'x']},
 {'ext': '_RRUU', 'def': ['UPPER', 'UNIT', 'CONJ'], 'undef': ['TRANSA'], 'for': ['c', 'z', 'x']},
 {'ext': '_RRUN', 'def': ['UPPER', 'CONJ'], 'undef': ['TRANSA', 'UNIT'], 'for': ['c', 'z', 'x']},
 {'ext': '_RRLU', 'def': ['UNIT', 'CONJ'], 'undef': ['TRANSA', 'UPPER'], 'for': ['c', 'z', 'x']},
 {'ext': '_RRLN', 'def': ['CONJ'], 'undef': ['TRANSA', 'UPPER', 'UNIT'], 'for': ['c', 'z', 'x']},
 {'ext': '_RCUU', 'def': ['TRANSA', 'UPPER', 'UNIT', 'CONJ'], 'undef': [], 'for': ['c', 'z', 'x']},
 {'ext': '_RCUN', 'def': ['TRANSA', 'UPPER', 'CONJ'], 'undef': ['UNIT'], 'for': ['c', 'z', 'x']},
 {'ext': '_RCLU', 'def': ['TRANSA', 'UNIT', 'CONJ'], 'undef': ['UPPER'], 'for': ['c', 'z', 'x']},
 {'ext': '_RCLN', 'def': ['TRANSA', 'CONJ'], 'undef': ['UPPER', 'UNIT'], 'for': ['c', 'z', 'x']},

 # symm
]

 # cc -c -O2 -DSMALL_MATRIX_OPT -DMAX_STACK_ALLOC=2048 -Wall -m64 -DF_INTERFACE_GFORT -fPIC -DSMP_SERVER -DNO_WARMUP -DMAX_CPU_NUMBER=12 -DMAX_PARALLEL_NUMBER=1 -DBUILD_SINGLE=1 -DBUILD_DOUBLE=1 -DBUILD_COMPLEX=1 -DBUILD_COMPLEX16=1 -DVERSION=\"0.3.26.dev\" -msse3 -mssse3 -msse4.1 -mavx -mavx2 -mavx2 -UASMNAME -UASMFNAME -UNAME -UCNAME -UCHAR_NAME -UCHAR_CNAME -DASMNAME=strmm_RTUU -DASMFNAME=strmm_RTUU_ -DNAME=strmm_RTUU_ -DCNAME=strmm_RTUU -DCHAR_NAME=\"strmm_RTUU_\" -DCHAR_CNAME=\"strmm_RTUU\" -DNO_AFFINITY -I../.. -UDOUBLE  -UCOMPLEX -UCOMPLEX -UDOUBLE -DTRANSA -DUPPER -DUNIT trmm_R.c -o strmm_RTUU.o

symb_defs = {
    '?amax':  {'def': ['USE_ABS'], 'undef': ['USE_MIN']},
    '?amin':  {'def': ['USE_ABS', 'USE_MIN']},
    'i?max':  {'undef': ['USE_ABS', 'USE_MIN']},
    'i?amax': {'def': ['USE_ABS'], 'undef': ['USE_MIN']},
    'i?amin': {'def': ['USE_ABS', 'USE_MIN']},
    'i?min':  {'def': ['USE_MIN'], 'undef': ['USE_ABS']},
    '?max':   {'undef': ['USE_ABS', 'USE_MIN']},
    '?min':   {'def': ['USE_MIN'], 'undef': ['USE_ABS']},
    '?axpyc': {'def': ['CONJ']},
    '?dotu':  {'undef': ['CONJ']},
    '?dotc':  {'def': ['CONJ']},
    '?geru':  {'undef': ['CONJ']},
    '?gerc':  {'def': ['CONJ']},
    '?hemm':  {'def': ['HEMM']},
    '?herk':  {'def': ['HEMM']},
    '?her2k':  {'def': ['HEMM']},
    '?gemm3m':  {'def': ['GEMM3M']},
    '?symm3m':  {'def': ['GEMM3M']},
    '?hemm3m':  {'def': ['HEMM', 'GEMM3M']},
    '?her_thread':  {'def': ['HER']},
    '?her2_thread':  {'def': ['HER']},
    '?hpr_thread':  {'def': ['HEMV']},
    '?bgemm':  {'def': ['HALF']},
    'cblas_?dotu_sub':  {'def': ['CBLAS', 'FORCE_USE_STACK'], 'undef': ['CONJ']},
    'cblas_?dotc_sub':  {'def': ['CBLAS', 'FORCE_USE_STACK', 'CONJ']},
}

# Ignoring other hostarch checks and conflicts for arch in BSD for now
_inc = [include_directories('.')]
# subdir('lapack-netlib')
subdir('interface')
subdir('driver/level2')
subdir('driver/level3')
subdir('kernel')

_openblas = static_library('openblas',
                    link_whole: [ _interface, _l2_driver, _l3_driver, _kern ])
