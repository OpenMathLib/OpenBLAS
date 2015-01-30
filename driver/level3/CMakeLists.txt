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
  # subtract 1 since we will iterate from 0 to num_combos
  math(EXPR num_combos "(${num_combos} << ${list_count}) - 1")
  set(LIST_OUT "")
  foreach (c RANGE 0 ${num_combos})
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
  list(APPEND LIST_OUT " ") # Empty set is a valic combination, but CMake isn't appending the empty string for some reason, use a space
  set(LIST_OUT ${LIST_OUT} PARENT_SCOPE)
endfunction ()

# generates object files for each of the sources for each of the combinations of the preprocessor definitions passed in
function(GenerateObjects sources_in defines_in)
  AllCombinations("${defines_in}")
  set(define_combos ${LIST_OUT})
  foreach (source_file ${sources_in})
    foreach (def_combo ${define_combos})

      # replace colon separated list with semicolons, this turns it into a CMake list that we can use foreach with
      string(REPLACE ":" ";" def_combo ${def_combo})

      # build a unique variable name for this obj file by picking two letters from the defines (can't use one in this case)
      set(obj_name "")
      foreach (combo_elem ${def_combo})
        string(REGEX MATCH "^[A-Z][A-Z]" letter ${combo_elem})
        set(obj_name "${obj_name}${letter}")
      endforeach ()

      # parse file name
      string(REGEX MATCH "^[a-zA-Z_]+" source_name ${source_file})
      string(TOUPPER ${source_name} source_name)

      # prepend the uppercased file name to the obj name
      set(obj_name "${source_name}_${obj_name}_OBJS")

      # now add the object and set the defines
      add_library(${obj_name} OBJECT ${source_file})
      if (NOT "${def_combo}" STREQUAL " ") # using space as the empty set
        set_target_properties(${obj_name} PROPERTIES COMPILE_DEFINITIONS "${def_combo}")
      endif ()
    endforeach ()
  endforeach ()
endfunction ()

# these sources are compiled with combinations of TRANS, UPPER, and UNIT, for 32 combinations total
set(TRM_SOURCES trmm_L.c trmm_R.c trsm_L.c trsm_R.c)
set(TRM_DEFINES TRANS UPPER UNIT)
GenerateObjects("${TRM_SOURCES}" "${TRM_DEFINES}")

# TODO: also need to set NN for all these objs (add param to GenerateObjects for defines that apply to all
GenerateObjects("symm_k.c" "LOWER;RSIDE")

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
