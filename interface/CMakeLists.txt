
include_directories(${CMAKE_SOURCE_DIR})

set(BLAS1_SOURCES
  axpy.c swap.c
  copy.c scal.c
  dot.c
  asum.c nrm2.c
  rot.c rotg.c rotm.c rotmg.c
  axpby.c
)

# TODO: USE_NETLIB_GEMV shoudl switch gemv.c to netlib/*gemv.f
set(BLAS2_SOURCES
  gemv.c ger.c
  trsv.c trmv.c symv.c
  syr.c syr2.c gbmv.c
  sbmv.c spmv.c
  spr.c spr2.c
  tbsv.c tbmv.c
  tpsv.c tpmv.c
)

set(BLAS3_SOURCES
  gemm.c symm.c
  trsm.c syrk.c syr2k.c
  omatcopy.c imatcopy.c
)

# generate the BLAS objs once with and once without cblas
set (CBLAS_FLAGS "")

if (NOT DEFINED NO_FBLAS)
  list(APPEND CBLAS_FLAGS 0)
endif ()

if (NOT DEFINED NO_CBLAS)
  list(APPEND CBLAS_FLAGS 1)
endif ()

foreach (CBLAS_FLAG ${CBLAS_FLAGS})

  GenerateNamedObjects("${BLAS1_SOURCES}" "" "" ${CBLAS_FLAG})
  GenerateNamedObjects("${BLAS2_SOURCES}" "" "" ${CBLAS_FLAG})
  GenerateNamedObjects("${BLAS3_SOURCES}" "" "" ${CBLAS_FLAG})

  # trmm is trsm with a compiler flag set
  GenerateNamedObjects("trsm.c" "TRMM" "trmm" ${CBLAS_FLAG})

  # max and imax are compiled 4 times
  GenerateNamedObjects("max.c" "" "" ${CBLAS_FLAG})
  GenerateNamedObjects("max.c" "USE_ABS" "amax" ${CBLAS_FLAG})
  GenerateNamedObjects("max.c" "USE_ABS;USE_MIN" "amin" ${CBLAS_FLAG})
  GenerateNamedObjects("max.c" "USE_MIN" "min" ${CBLAS_FLAG})

  GenerateNamedObjects("imax.c" "" "i*max" ${CBLAS_FLAG})
  GenerateNamedObjects("imax.c" "USE_ABS" "i*amax" ${CBLAS_FLAG})
  GenerateNamedObjects("imax.c" "USE_ABS;USE_MIN" "i*amin" ${CBLAS_FLAG})
  GenerateNamedObjects("imax.c" "USE_MIN" "i*min" ${CBLAS_FLAG})

endforeach ()

if (NOT DEFINED NO_LAPACK)
  set(LAPACK_SOURCES
    lapack/getrf.c lapack/getrs.c lapack/potrf.c lapack/getf2.c
    lapack/potf2.c lapack/laswp.c lapack/gesv.c lapack/lauu2.c
    lapack/lauum.c lapack/trti2.c lapack/trtri.c
  )
  GenerateNamedObjects("${LAPACK_SOURCES}")
endif ()

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

