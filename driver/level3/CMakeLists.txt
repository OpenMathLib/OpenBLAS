include_directories(${CMAKE_SOURCE_DIR})

set(USE_GEMM3M 0)

if (DEFINED ARCH)
  if (${ARCH} STREQUAL "x86")
    set(USE_GEMM3M 1)
  endif ()

  if (${ARCH} STREQUAL "x86_64")
    set(USE_GEMM3M 1)
  endif ()

  if (${ARCH} STREQUAL "ia64")
    set(USE_GEMM3M 1)
  endif ()

  if (${ARCH} STREQUAL "MIPS")
    set(USE_GEMM3M 1)
  endif ()
endif ()

# N.B. In the original makefile there was a BLOCKS define used in the compilation of these files but I don't see any evidence of it being set anywhere. -hpa

# loop through gemm.c defines
set(GEMM_DEFINES NN NT TN TT)
foreach (GEMM_DEFINE ${GEMM_DEFINES})
  string(TOLOWER ${GEMM_DEFINE} GEMM_DEFINE_LC)
  GenerateNamedObjects("gemm.c" "DOUBLE" "${GEMM_DEFINE}" "gemm_${GEMM_DEFINE_LC}" 0)
  if (SMP AND NOT USE_SIMPLE_THREADED_LEVEL3)
    GenerateNamedObjects("gemm.c" "DOUBLE" "${GEMM_DEFINE};THREADED_LEVEL3" "gemm_thread_${GEMM_DEFINE_LC}" 0)
  endif ()
endforeach ()

GenerateCombinationObjects("trmm_L.c;trmm_R.c;trsm_L.c;trsm_R.c" "TRANS;UPPER;UNIT" "N;L;N" "DOUBLE" "" 0)
GenerateCombinationObjects("symm_k.c" "RSIDE;LOWER" "L;U" "DOUBLE" "NN" 1)
GenerateCombinationObjects("syrk_k.c;syr2k_k.c" "LOWER;TRANS" "U;N" "DOUBLE" "" 1)
GenerateCombinationObjects("syrk_kernel.c;syr2k_kernel.c" "LOWER" "U" "DOUBLE" "" 2)

if (SMP)

  # N.B. these do NOT have a float type (e.g. DOUBLE) defined!
  GenerateNamedObjects("gemm_thread_m.c;gemm_thread_n.c;gemm_thread_mn.c;gemm_thread_variable.c;syrk_thread.c" "" "" "" 0)

  if (NOT USE_SIMPLE_THREADED_LEVEL3)
    GenerateCombinationObjects("syrk_k.c" "LOWER;TRANS" "U;N" "DOUBLE" "THREADED_LEVEL3" 2 "syrk_thread")
    GenerateCombinationObjects("symm_k.c" "RSIDE;LOWER" "L;U" "DOUBLE" "THREADED_LEVEL3;NN" 2 "symm_thread")
  endif ()
endif ()

#HPLOBJS =
#	dgemm_nn.c dgemm_nt.c dgemm_tn.c dgemm_tt.c
#	dtrsm_LNUU.c dtrsm_LNUN.c dtrsm_LNLU.c dtrsm_LNLN.c
#	dtrsm_LTUU.c dtrsm_LTUN.c dtrsm_LTLU.c dtrsm_LTLN.c
#	dtrsm_RNUU.c dtrsm_RNUN.c dtrsm_RNLU.c dtrsm_RNLN.c
#	dtrsm_RTUU.c dtrsm_RTUN.c dtrsm_RTLU.c dtrsm_RTLN.c
#
#if (USE_SIMPLE_THREADED_LEVEL3)
#  HPLOBJS += dgemm_thread_nn.c dgemm_thread_nt.c
#	   dgemm_thread_tn.c dgemm_thread_tt.c
#endif
#

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

