
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

# first compile all the objects that don't need specific preprocessor defines
add_library(DBLAS_NONE OBJECT
  gbmv_k.c # gbmv_N
  ${UL_SOURCES}
  ${NU_SOURCES}
)

# then do objects with transpose/triangular/etc definitions

# objects that need TRANS set
add_library(DBLAS_T OBJECT gbmv_k.c ${NU_SOURCES})
set_target_properties(DBLAS_T PROPERTIES COMPILE_DEFINITIONS TRANS)

# objects that need LOWER set
add_library(DBLAS_L OBJECT ${UL_SOURCES})
set_target_properties(DBLAS_L PROPERTIES COMPILE_DEFINITIONS LOWER)

# objects that need UNIT set
add_library(DBLAS_U OBJECT ${NU_SOURCES})
set_target_properties(DBLAS_U PROPERTIES COMPILE_DEFINITIONS UNIT)

# objects that need TRANS and UNIT set
add_library(DBLAS_TU OBJECT ${NU_SOURCES})
set_target_properties(DBLAS_TU PROPERTIES COMPILE_DEFINITIONS UNIT)
set_target_properties(DBLAS_TU PROPERTIES COMPILE_DEFINITIONS TRANS)

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

set(DBLAS_TARGETS DBLAS_NONE DBLAS_T DBLAS_L DBLAS_U DBLAS_TU)

foreach (${TARGET} ${DBLAS_TARGETS})
  set_target_properties(${TARGET} PROPERTIES COMPILE_DEFINITIONS DOUBLE)
endforeach ()

