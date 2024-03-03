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

# Makefile.system
# Ignoring all the hostarch checks and conflits for arch in BSD for now
_inc = include_directories('.')
# subdir('lapack-netlib')
subdir('interface')

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

if host_machine.system() == 'linux'
    no_affinity = true
else
    no_affinity = false
endif

# Example for handling options:
if build_without_lapack
    # configure build to exclude LAPACK and LAPACKE
endif
