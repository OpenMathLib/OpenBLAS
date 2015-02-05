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
  GenerateNamedObjects("gemm.c" "DOUBLE" "${GEMM_DEFINE}" "gemm_${GEMM_DEFINE_LC}" 0 "" "")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
endforeach ()

AllCombinations("TRANS;UPPER;UNIT")
set(define_combos ${LIST_OUT})
foreach (def_combo ${define_combos})
  # replace colon separated list with semicolons, this turns it into a CMake list that we can use foreach with
  string(REPLACE ":" ";" def_combo ${def_combo})

  # build a unique variable name for this obj file by picking two letters from the defines (can't use one in this case)
  set(obj_name "")
  foreach (combo_elem ${def_combo})
    string(REGEX MATCH "^[A-Z][A-Z]" letter ${combo_elem})
    set(obj_name "${obj_name}${letter}")
  endforeach ()
endforeach ()

GenerateCombinationObjects("trmm_L.c;trmm_R.c;trsm_L.c;trsm_R.c" "TRANS;UPPER;UNIT" "DOUBLE" "" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateCombinationObjects("symm_k.c" "LOWER;RSIDE" "DOUBLE" "NN" 1)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateCombinationObjects("syrk_k.c;syr2k_k.c" "LOWER;TRANS" "DOUBLE" "" 1)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateCombinationObjects("syrk_kernel.c;syr2k_kernel.c" "LOWER" "DOUBLE" "" 1)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

#if (SMP)
#
#  COMMONOBJS  += gemm_thread_m.c gemm_thread_n.c gemm_thread_mn.c gemm_thread_variable.c
#  COMMONOBJS  += syrk_thread.c
#
#  if (USE_SIMPLE_THREADED_LEVEL3)
#    DBLASOBJS    += dgemm_thread_nn.c dgemm_thread_nt.c dgemm_thread_tn.c dgemm_thread_tt.c
#    DBLASOBJS    += dsymm_thread_LU.c dsymm_thread_LL.c dsymm_thread_RU.c dsymm_thread_RL.c
#    DBLASOBJS    += dsyrk_thread_UN.c dsyrk_thread_UT.c dsyrk_thread_LN.c dsyrk_thread_LT.c
#
#  endif ()
#endif ()
#
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

