# Relevant groups and variables:
#
# _install_aux -- Auxiliary function library
# _int_cpu_time -- Timing library, Fortran only
#

prec = get_option('realkind')
build_single = get_option('build_single')
build_double = get_option('build_double')
build_complex = get_option('build_complex')
build_complex16 = get_option('build_complex16')
build_all_prec = get_option('build_all_prec')
use_c_lapack = get_option('use_c_lapack')

# C / Fortran
_install_aux_srcs = []

if not use_c_lapack
  _int_cpu_srcs = [] # Fortran only
  _install_aux_srcs += ['ilaver.f']

  if prec == 'd' or build_double or build_all_prec
    _int_cpu_srcs += [
      'dsecnd_INT_CPU_TIME.f',
    ]
    _install_aux_srcs += [
      'dlamch.f',
      'droundup_lwork.f',
    ]
  endif

  if prec == 's' or build_single or build_all_prec
    _int_cpu_srcs += [
      'second_INT_CPU_TIME.f',
    ]
    _install_aux_srcs += [
      'slamch.f',
      'sroundup_lwork.f',
    ]
  endif
  _int_cpu = library('_int_cpu',
                     sources: _int_cpu_srcs,
                     version: lapack_version,
                     soversion: lapack_major_version,
                     install: true)
else
  _install_aux_srcs += ['ilaver.c']

  if prec == 'd' or build_double or build_all_prec
    _install_aux_srcs += [
      'dlamch.c'
    ]
  endif

  if prec == 's' or build_single or build_all_prec
    _install_aux_srcs += [
      'slamch.c'
    ]
  endif

endif

_install_aux = library('_install_aux',
                       sources: install_aux_srcs,
                       version: lapack_version,
                       soversion: lapack_major_version,
                       install: true)
