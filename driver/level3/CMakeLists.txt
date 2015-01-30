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
  add_library(GEMM_${GEMM_DEFINE}_OBJS OBJECT gemm.c)
  set_target_properties(GEMM_${GEMM_DEFINE}_OBJS PROPERTIES COMPILE_DEFINITIONS ${GEMM_DEFINE})
endforeach ()

# Returns all combinations of the input list, as a list with colon-separated combinations
# E.g. input of A B C returns A B C A:B A:C B:C
# N.B. The input is meant to be a list, and to past a list to a function in CMake you must quote it (e.g. AllCombinations("${LIST_VAR}")).
function(AllCombinations list_in)
  list(LENGTH list_in list_count)
  set(num_combos 1)
  math(EXPR num_combos "${num_combos} << ${list_count}")
  set(LIST_OUT "")
  foreach (c RANGE ${num_combos})
    set(current_combo "")
    # this is a little ridiculous just to iterate through a list w/ indices
    math(EXPR last_list_index "${list_count} - 1")
    foreach (list_index RANGE 0 ${last_list_index})
      math(EXPR bit "1 << ${list_index}")
      math(EXPR combo_has_bit "${c} & ${bit}")
      list(GET list_in ${list_index} list_elem)
      if (combo_has_bit)
        if (current_combo)
          set(current_combo "${current_combo}:${list_elem}")
        else ()
          set(current_combo ${list_elem})
        endif ()
      endif ()
    endforeach ()
    list(APPEND LIST_OUT ${current_combo})
  endforeach ()
  set(LIST_OUT ${LIST_OUT} PARENT_SCOPE)
endfunction ()

# these sources are compiled with combinations of TRANS, UPPER, and UNIT, for 32 combinations total
set(TRM_SOURCES trmm_L.c trmm_R.c trsm_L.c trsm_R.c)
AllCombinations("TRANS UPPER UNIT")
set(TRM_DEFINE_COMBOS LIST_OUT)
message(STATUS "alcombos result: ${LIST_OUT}")
foreach (TRM_SOURCE ${TRM_SOURCES})
  foreach (TRM_DEFINES ${TRM_DEFINE_COMBOS})
    string(REGEX MATCH "[a-z]+_[LR]" TRM_NAME ${TRM_SOURCE})
    string(TOUPPER ${TRM_NAME} TRM_NAME)
    # TODO: TRM_DEFINES is a colon-separated list of defines to set for this object - need to parse it and set them using set_target_properties, and also come up with a unique id for the lib name (e.g. first letter of each define, so TRANS UPPER UNIT is TUU)
    #add_library(${TRM_NAME}_${TRM_DEFINE}_OBJS OBJECT ${TRM_SOURCE})
    #set_target_properties(${TRM_NAME}_${TRM_DEFINE}_OBJS PROPERTIES COMPILE_DEFINITIONS ${TRM_DEFINE})
  endforeach ()
endforeach ()

#	dsymm_LU.c dsymm_LL.c dsymm_RU.c dsymm_RL.c
#	dsyrk_UN.c dsyrk_UT.c dsyrk_LN.c dsyrk_LT.c
#	dsyr2k_UN.c dsyr2k_UT.c dsyr2k_LN.c dsyr2k_LT.c
#	dsyrk_kernel_U.c  dsyrk_kernel_L.c
#	dsyr2k_kernel_U.c dsyr2k_kernel_L.c

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
