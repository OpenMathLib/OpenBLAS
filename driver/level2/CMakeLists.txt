
include_directories(${CMAKE_SOURCE_DIR})

# sources that need to be compiled twice, once with no flags and once with LOWER
set(UL_SOURCES
  sbmv_k.c
  spmv_k.c
  spr_k.c
  spr2_k.c
  syr_k.c
  syr2_k.c
)

# sources that need to be compiled several times, for UNIT, TRANS
set(NU_SOURCES
  tbmv_U.c
  tbsv_U.c
  tpmv_U.c
  tpsv_U.c
  trmv_U.c
  trsv_U.c
  tbmv_L.c
  tbsv_L.c
  tpmv_L.c
  tpsv_L.c
  trmv_L.c
  trsv_L.c
)

# objects that need LOWER set
GenerateCombinationObjects("${UL_SOURCES}" "LOWER" "U" "DOUBLE" "" 1)
list(APPEND DBLAS_OBJS ${COMBO_OBJ_LIST_OUT})

# objects that need TRANS and UNIT set
# N.B. BLAS wants to put the U/L from the filename in the *MIDDLE* because of course why not have a different naming scheme for every single object -hpa
GenerateCombinationObjects("${NU_SOURCES}" "TRANS;UNIT" "N;N" "DOUBLE" "" 3)
list(APPEND DBLAS_OBJS ${COMBO_OBJ_LIST_OUT})

# gbmv uses a lowercase n and t. WHY? WHO KNOWS!
GenerateNamedObjects("gbmv_k.c" "DOUBLE" "" "gbmv_n" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("gbmv_k.c" "DOUBLE" "TRANS" "gbmv_t" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

#if (DEFINED SMP)
#  add_library(DBLASOBJS_SMP
#    dgemv_thread_n.c	dgemv_thread_t.c
#    dger_thread.c
#    dsymv_thread_U.c	dsymv_thread_L.c
#    dsyr_thread_U.c		dsyr_thread_L.c
#    dsyr2_thread_U.c	dsyr2_thread_L.c
#    dspr_thread_U.c		dspr_thread_L.c
#    dspr2_thread_U.c	dspr2_thread_L.c
#    dtrmv_thread_NUU.c	dtrmv_thread_NUN.c
#    dtrmv_thread_NLU.c	dtrmv_thread_NLN.c
#    dtrmv_thread_TUU.c	dtrmv_thread_TUN.c
#    dtrmv_thread_TLU.c	dtrmv_thread_TLN.c
#    dspmv_thread_U.c	dspmv_thread_L.c
#    dtpmv_thread_NUU.c	dtpmv_thread_NUN.c
#    dtpmv_thread_NLU.c	dtpmv_thread_NLN.c
#    dtpmv_thread_TUU.c	dtpmv_thread_TUN.c
#    dtpmv_thread_TLU.c	dtpmv_thread_TLN.c
#    dgbmv_thread_n.c	dgbmv_thread_t.c
#    dsbmv_thread_U.c	dsbmv_thread_L.c
#    dtbmv_thread_NUU.c	dtbmv_thread_NUN.c
#    dtbmv_thread_NLU.c	dtbmv_thread_NLN.c
#    dtbmv_thread_TUU.c	dtbmv_thread_TUN.c
#    dtbmv_thread_TLU.c	dtbmv_thread_TLN.c
#  )
#endif ()

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

