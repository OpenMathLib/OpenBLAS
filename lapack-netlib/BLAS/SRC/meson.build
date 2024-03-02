# Derived from the CMakeLists.txt
#
# Relevant groups and variables:
#
# _allblas -- Auxiliary routines for Level 2 and 3 BLAS
#
# Level 1 BLAS #
#
# _dblas1 -- Double precision real BLAS 1 routines
# _zblas1 -- Double precision complex BLAS 1 routines
# _zb1aux -- D.P. real BLAS 1 routines called by d.p. complex routines
#
# _sblas1 -- Single precision real BLAS routines
# _cblas1 -- Single precision complex BLAS routines
# _cb1aux -- Real BLAS routines called by complex routines
#
# Level 2 BLAS #
#
# _dblas2 -- Double precision real BLAS 2 routines
# _zblas2 -- Double precision complex BLAS 2 routines
#
# _sblas2 -- Single precision real BLAS2 routines
# _cblas2 -- Single precision complex BLAS2 routines
#
# Level 3 BLAS #
#
# _dblas3 -- Double precision real BLAS 3 routines
# _zblas3 -- Double precision complex BLAS 3 routines
#
# _sblas3 -- Single precision real BLAS3 routines
# _cblas3 -- Single precision complex BLAS3 routines

# _allblas -- Auxiliary routines for Level 2 and 3 BLAS
_allblas = library('_allblas',
                   sources: [ 'lsame.f', 'xerbla.f', 'xerbla_array.f' ])


# All other sources
_blas_netlib_srcs = []

# Level 1 BLAS
# _dblas1 -- Double precision real BLAS 1 routines
_dblas1 = [
  'idamax.f',
  'dasum.f',
  'daxpy.f',
  'dcopy.f',
  'ddot.f',
  'dnrm2.f',
  'drot.f',
  'drotg.f',
  'dscal.f',
  'dsdot.f',
  'dswap.f',
  'drotmg.f',
  'drotm.f',
]

# _zblas1 -- Double precision complex BLAS 1 routines
_zblas1 = [
  'dcabs1.f',
  'dzasum.f',
  'dznrm2.f',
  'izamax.f',
  'zaxpy.f',
  'zcopy.f',
  'zdotc.f',
  'zdotu.f',
  'zdscal.f',
  'zrotg.f',
  'zscal.f',
  'zswap.f',
  'zdrot.f',
]

# _zb1aux -- D.P. real BLAS routines called by d.p. complex routines
_zb1aux = [
  'idamax.f',
  'dasum.f',
  'daxpy.f',
  'dcopy.f',
  'dnrm2.f',
  'dscal.f',
]

# _sblas1 -- Single precision real BLAS routines
_sblas1 = [
  'isamax.f',
  'sasum.f',
  'saxpy.f',
  'scopy.f',
  'sdot.f',
  'snrm2.f',
  'srot.f',
  'srotg.f',
  'sscal.f',
  'sswap.f',
  'sdsdot.f',
  'srotmg.f',
  'srotm.f',
]

# _cblas1 -- Single precision complex BLAS routines
_cblas1 = [
  'scabs1.f',
  'scasum.f',
  'scnrm2.f',
  'icamax.f',
  'caxpy.f',
  'ccopy.f',
  'cdotc.f',
  'cdotu.f',
  'csscal.f',
  'crotg.f',
  'cscal.f',
  'cswap.f',
  'csrot.f',
]

# _cb1aux -- Real BLAS routines called by complex routines
_cb1aux = [
  'isamax.f',
  'sasum.f',
  'saxpy.f',
  'scopy.f',
  'snrm2.f',
  'sscal.f'
]

# Level 2 BLAS
# _dblas2 -- Double precision real BLAS2 routines
_dblas2 = [
  'dgemv.f',
  'dgbmv.f',
  'dsymv.f',
  'dsbmv.f',
  'dspmv.f',
  'dtrmv.f',
  'dtbmv.f',
  'dtpmv.f',
  'dtrsv.f',
  'dtbsv.f',
  'dtpsv.f',
  'dger.f',
  'dsyr.f',
  'dspr.f',
  'dsyr2.f',
  'dspr2.f',
]

# _zblas2 -- Double precision complex BLAS2 routines
_zblas2 = [
  'zgemv.f',
  'zgbmv.f',
  'zhemv.f',
  'zhbmv.f',
  'zhpmv.f',
  'ztrmv.f',
  'ztbmv.f',
  'ztpmv.f',
  'ztrsv.f',
  'ztbsv.f',
  'ztpsv.f',
  'zgerc.f',
  'zgeru.f',
  'zher.f',
  'zhpr.f',
  'zher2.f',
  'zhpr2.f',
]

# _sblas2 -- Single precision real BLAS2 routines
_sblas2 = [
  'sgemv.f',
  'sgbmv.f',
  'ssymv.f',
  'ssbmv.f',
  'sspmv.f',
  'strmv.f',
  'stbmv.f',
  'stpmv.f',
  'strsv.f',
  'stbsv.f',
  'stpsv.f',
  'sger.f',
  'ssyr.f',
  'sspr.f',
  'ssyr2.f',
  'sspr2.f',
]

# _cblas2 -- Single precision complex BLAS2 routines
_cblas2 = [
  'cgemv.f',
  'cgbmv.f',
  'chemv.f',
  'chbmv.f',
  'chpmv.f',
  'ctrmv.f',
  'ctbmv.f',
  'ctpmv.f',
  'ctrsv.f',
  'ctbsv.f',
  'ctpsv.f',
  'cgerc.f',
  'cgeru.f',
  'cher.f',
  'chpr.f',
  'cher2.f',
  'chpr2.f',
]

# Level 3 BLAS

# _dblas3 -- Double precision real BLAS3 routines
_dblas3 = [
  'dgemm.f',
  'dsymm.f',
  'dsyrk.f',
  'dsyr2k.f',
  'dtrmm.f',
  'dtrsm.f',
]

# _zblas3 -- Double precision complex BLAS3 routines
_zblas3 = [
  'zgemm.f',
  'zsymm.f',
  'zsyrk.f',
  'zsyr2k.f',
  'ztrmm.f',
  'ztrsm.f',
  'zhemm.f',
  'zherk.f',
  'zher2k.f',
]

# _sblas3 -- Single precision real BLAS3 routines
_sblas3 =[
  'sgemm.f',
  'ssymm.f',
  'ssyrk.f',
  'ssyr2k.f',
  'strmm.f',
  'strsm.f',
]

# _cblas3 -- Single precision complex BLAS3 routines
_cblas3 = [
  'cgemm.f',
  'csymm.f',
  'csyrk.f',
  'csyr2k.f',
  'ctrmm.f',
  'ctrsm.f',
  'chemm.f',
  'cherk.f',
  'cher2k.f',
]

# Start making the blas target
if build_complex or build_all_prec
  _blas_netlib_srcs += _cblas1 + _cb1aux + _cblas2 + _cblas3
endif
if build_complex16 or build_all_prec
  _blas_netlib_srcs += _zblas1 + _zb1aux + _zblas2 + _zblas3
endif

if prec == 'd' or build_double or build_all_prec
  _blas_netlib_srcs += _dblas1 + _dblas2 + _dblas3
endif

if prec == 's' or build_single or build_all_prec
  _blas_netlib_srcs += _sblas1 + _sblas2 + _sblas3
endif

# Create the blas library
netlib_blas = library(_netlib_blas_name,
               sources: _blas_netlib_srcs,
               link_with: _allblas,
               version: lapack_version,
               soversion: lapack_major_version,
               install: true)
