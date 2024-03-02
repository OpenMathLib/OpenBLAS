# Relevant groups and variables:
#
# _install_aux -- Auxiliary function library
# _int_cpu_time -- Timing library, Fortran only
#
# C / Fortran
_install_aux_srcs = []
_ia_linkwith = []

if not use_c_lapack
  _int_cpu_srcs = [] # Fortran only
  _install_aux_srcs += ['ilaver.f']
  _ia_linkwith = _allblas # Includes lsame.f
  # TODO: Delete the one in install in favor of the one from BLAS
  # The meson build uses the BLAS variant anyway
  # lsame.f is the same modulo formatting changes as the one in blas/src/lsame.f
  # lsame.f claims to be from lapack 3.7.0 and has more whitespace but is
  # otherwise identical to the reference from blas 3.1

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
  _install_aux_srcs += ['ilaver.c', 'lsame.c']

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
                       sources: _install_aux_srcs,
                       link_with: _ia_linkwith,
                       version: lapack_version,
                       soversion: lapack_major_version,
                       install: true)
