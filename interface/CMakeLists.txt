
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


if (NOT DEFINED NO_FBLAS)

  GenerateNamedObjects("${BLAS1_SOURCES}" "DOUBLE" "" "")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("${BLAS2_SOURCES}" "DOUBLE" "" "")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("${BLAS3_SOURCES}" "DOUBLE" "" "")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

  # trmm is trsm with a compiler flag set
  GenerateNamedObjects("trsm.c" "DOUBLE" "TRMM" "trmm")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

  # max and imax are compiled 4 times
  GenerateNamedObjects("max.c" "DOUBLE" "" "")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("max.c" "DOUBLE" "USE_ABS" "amax")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("max.c" "DOUBLE" "USE_ABS;USE_MIN" "amin")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("max.c" "DOUBLE" "USE_MIN" "min")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

  GenerateNamedObjects("imax.c" "DOUBLE" "" "i*max")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("imax.c" "DOUBLE" "USE_ABS" "i*amax")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("imax.c" "DOUBLE" "USE_ABS;USE_MIN" "i*amin")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
  GenerateNamedObjects("imax.c" "DOUBLE" "USE_MIN" "i*min")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

endif ()

if (NOT DEFINED NO_CBLAS)

  add_library(ISAMAX_OBJ OBJECT imax.c)
  set_target_properties(ISAMAX_OBJ PROPERTIES COMPILE_DEFINITIONS "CBLAS;USE_ABS")

  add_library(CDBLAS1_OBJS OBJECT ${BLAS1_SOURCES})
  add_library(CDBLAS2_OBJS OBJECT ${BLAS2_SOURCES})
  add_library(CDBLAS3_OBJS OBJECT ${BLAS3_SOURCES})

  # trmm is trsm with a compiler flag set
  add_library(CTRMM_OBJ OBJECT trsm.c)
  set_target_properties(CTRMM_OBJ PROPERTIES COMPILE_DEFINITIONS "CBLAS;TRMM")

  set_target_properties(CDBLAS1_OBJS PROPERTIES COMPILE_DEFINITIONS "CBLAS")
  set_target_properties(CDBLAS2_OBJS PROPERTIES COMPILE_DEFINITIONS "CBLAS")
  set_target_properties(CDBLAS3_OBJS PROPERTIES COMPILE_DEFINITIONS "CBLAS")

  list(APPEND DBLAS_OBJS "CDBLAS1_OBJS;CDBLAS2_OBJS;CDBLAS3_OBJS;ISAMAX_OBJ;CTRMM_OBJ")
endif ()

if (NOT DEFINED NO_LAPACK)
  add_library(DLAPACK_OBJS OBJECT
    lapack/getrf.c lapack/getrs.c lapack/potrf.c lapack/getf2.c
    lapack/potf2.c lapack/laswp.c lapack/gesv.c lapack/lauu2.c
    lapack/lauum.c lapack/trti2.c lapack/trtri.c
  )
  list(APPEND DBLAS_OBJS "DLAPACK_OBJS")
endif ()

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

