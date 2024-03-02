prec = get_option('realkind')
if prec == 'd'
  if not get_option('use_c_lapack')
    LINSTALL = files(
      'droundup_lwork.f',
      'lsame.f',
    )
  else
    LINSTALL = files(
      'dlamch.c',
      'lsame.c',
    )
  endif
elif prec == 's'
  if not get_option('use_c_lapack')
    LINSTALL = files(
      'slamch.f',
      'sroundup_lwork.f',
      'lsame.f',
    )
  else
    LINSTALL = files(
      'slamch.c',
      'lsame.c',
    )
  endif
endif
