
include_directories(${CMAKE_SOURCE_DIR})

# TODO: Need to generate object files for S, D, C, Q and X - start with D for now.
#       The sources are the same, but there are additional preprocessor definitions depending on the precision (see Makefile.tail).

add_library(DBLAS1OBJS OBJECT
  axpy.c swap.c
  copy.c scal.c
  dot.c
  asum.c nrm2.c
  max.c # amax/min/amin compiled later from same source
  rot.c rotg.c rotm.c rotmg.c
  axpby.c
)

# N.B. The original Makefile passed in -UUSE_MIN and -UUSE_ABS (where appropriate), no way to do that at a source-level in cmake. REMOVE_DEFINITIONS removes a definition for the rest of the compilation.
add_library(AMAX_OBJ OBJECT max.c)
set_target_properties(AMAX_OBJ PROPERTIES COMPILE_DEFINITIONS "USE_ABS")
add_library(AMIN_OBJ OBJECT max.c)
set_target_properties(AMIN_OBJ PROPERTIES COMPILE_DEFINITIONS "USE_ABS;USE_MIN")
add_library(MIN_OBJ OBJECT max.c)
set_target_properties(MIN_OBJ PROPERTIES COMPILE_DEFINITIONS "USE_MIN")

# TODO: USE_NETLIB_GEMV shoudl switch gemv.c to netlib/*gemv.f
add_library(DBLAS2OBJS OBJECT
  gemv.c ger.c
  trsv.c trmv.c symv.c
  syr.c  syr2.c gbmv.c
  sbmv.c spmv.c
  spr.c  spr2.c
  tbsv.c tbmv.c
  tpsv.c tpmv.c
)

add_library(DBLAS3OBJS OBJECT
  gemm.c symm.c
  trsm.c syrk.c syr2k.c
  omatcopy.c imatcopy.c
)

# trmm is trsm with a compiler flag set
add_library(TRMM_OBJ OBJECT trsm.c)
set_target_properties(TRMM_OBJ PROPERTIES COMPILE_DEFINITIONS TRMM)

