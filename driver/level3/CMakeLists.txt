include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

# N.B. In the original makefile there was a BLOCKS define used in the compilation of these files but I don't see any evidence of it being set anywhere. -hpa

# loop through gemm.c defines
set(GEMM_DEFINES NN NT TN TT)
set(GEMM_COMPLEX_DEFINES RN CN RT CT NR TR RR CR NC TC RC CC)
foreach (GEMM_DEFINE ${GEMM_DEFINES})
  string(TOLOWER ${GEMM_DEFINE} GEMM_DEFINE_LC)
  GenerateNamedObjects("gemm.c" "${GEMM_DEFINE}" "gemm_${GEMM_DEFINE_LC}" 0)
  if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
    GenerateNamedObjects("gemm.c" "${GEMM_DEFINE};THREADED_LEVEL3" "gemm_thread_${GEMM_DEFINE_LC}" 0)
  endif ()
  if (BUILD_BFLOAT16)
    GenerateNamedObjects("gemm.c" "${GEMM_DEFINE}" "gemm_${GEMM_DEFINE_LC}" 0 "" "" false "BFLOAT16")
    GenerateNamedObjects("gemm.c" "${GEMM_DEFINE};BGEMM" "gemm_${GEMM_DEFINE_LC}" 0 "" "" false "BFLOAT16")
    if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
      GenerateNamedObjects("gemm.c" "${GEMM_DEFINE};THREADED_LEVEL3" "gemm_thread_${GEMM_DEFINE_LC}" 0 "" "" false "BFLOAT16")
      GenerateNamedObjects("gemm.c" "${GEMM_DEFINE};THREADED_LEVEL3;BGEMM" "gemm_thread_${GEMM_DEFINE_LC}" 0 "" "" false "BFLOAT16")
    endif ()
  endif ()
  if (BUILD_HFLOAT16)
    GenerateNamedObjects("gemm.c" "${GEMM_DEFINE}" "gemm_${GEMM_DEFINE_LC}" 0 "" "" false "HFLOAT16")
    if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
      GenerateNamedObjects("gemm.c" "${GEMM_DEFINE};THREADED_LEVEL3" "gemm_thread_${GEMM_DEFINE_LC}" 0 "" "" false "HFLOAT16")
    endif ()
  endif ()
endforeach ()

if ( BUILD_COMPLEX16 AND NOT  BUILD_DOUBLE)
foreach (GEMM_DEFINE ${GEMM_DEFINES})
  string(TOLOWER ${GEMM_DEFINE} GEMM_DEFINE_LC)
  GenerateNamedObjects("gemm.c" "${GEMM_DEFINE}" "gemm_${GEMM_DEFINE_LC}" 0 "" "" false "DOUBLE")
  if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
	  GenerateNamedObjects("gemm.c" "${GEMM_DEFINE};THREADED_LEVEL3" "gemm_thread_${GEMM_DEFINE_LC}" 0 "" "" false "DOUBLE")
  endif()
endforeach()
endif()
if ( BUILD_COMPLEX AND NOT  BUILD_SINGLE)
foreach (GEMM_DEFINE ${GEMM_DEFINES})
  string(TOLOWER ${GEMM_DEFINE} GEMM_DEFINE_LC)
  GenerateNamedObjects("gemm.c" "${GEMM_DEFINE}" "gemm_${GEMM_DEFINE_LC}" 0 "" "" false "SINGLE")
  if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
	  GenerateNamedObjects("gemm.c" "${GEMM_DEFINE};THREADED_LEVEL3" "gemm_thread_${GEMM_DEFINE_LC}" 0 "" "" false "SINGLE")
  endif()
endforeach()
endif()

set(TRMM_TRSM_SOURCES
  trmm_L.c
  trmm_R.c
  trsm_L.c
  trsm_R.c)

foreach(trmm_trsm_source ${TRMM_TRSM_SOURCES})
  string(REGEX MATCH "[a-z]+_[A-Z]+" op_name ${trmm_trsm_source})
  GenerateCombinationObjects("${trmm_trsm_source}" "UPPER;UNIT" "L;N" "" 0 "${op_name}N")
  GenerateCombinationObjects("${trmm_trsm_source}" "UPPER;UNIT" "L;N" "TRANSA" 0 "${op_name}T")
endforeach()

GenerateCombinationObjects("symm_k.c" "RSIDE;LOWER" "L;U" "NN" 1)
GenerateCombinationObjects("syrk_k.c" "LOWER;TRANS" "U;N" "" 1)
GenerateCombinationObjects("syr2k_k.c" "LOWER;TRANS" "U;N" "" 1)
GenerateCombinationObjects("syrk_kernel.c" "LOWER" "U" "" 2)

# GEMMT / GEMMTR: one driver per uplo x op(A) x op(B) (the same set of
# op combinations as GEMM), plus the conjugating variants of the
# triangular kernel, which only complex GEMMT needs.
set(GEMMT_REAL_TRANS N T)
set(GEMMT_COMPLEX_TRANS N T R C)
foreach (GEMMT_UPLO U L)
  set(GEMMT_UPLO_DEFS "")
  if (${GEMMT_UPLO} STREQUAL "L")
    set(GEMMT_UPLO_DEFS "LOWER")
  endif ()
  foreach (GEMMT_TA ${GEMMT_REAL_TRANS})
    foreach (GEMMT_TB ${GEMMT_REAL_TRANS})
      GenerateNamedObjects("gemmt_k.c" "${GEMMT_TA}${GEMMT_TB};${GEMMT_UPLO_DEFS}" "gemmt_${GEMMT_UPLO}${GEMMT_TA}${GEMMT_TB}" 0 "" "" false 1)
    endforeach ()
  endforeach ()
  foreach (GEMMT_TA ${GEMMT_COMPLEX_TRANS})
    foreach (GEMMT_TB ${GEMMT_COMPLEX_TRANS})
      GenerateNamedObjects("gemmt_k.c" "${GEMMT_TA}${GEMMT_TB};${GEMMT_UPLO_DEFS}" "gemmt_${GEMMT_UPLO}${GEMMT_TA}${GEMMT_TB}" 0 "" "" false 2)
    endforeach ()
  endforeach ()
  GenerateNamedObjects("syrk_kernel.c" "CONJA;${GEMMT_UPLO_DEFS}" "gemmt_kernel_${GEMMT_UPLO}CN" 0 "" "" false 2)
  GenerateNamedObjects("syrk_kernel.c" "CONJB;${GEMMT_UPLO_DEFS}" "gemmt_kernel_${GEMMT_UPLO}NC" 0 "" "" false 2)
  GenerateNamedObjects("syrk_kernel.c" "CONJA;CONJB;${GEMMT_UPLO_DEFS}" "gemmt_kernel_${GEMMT_UPLO}CC" 0 "" "" false 2)
endforeach ()
GenerateCombinationObjects("syr2k_kernel.c" "LOWER" "U" "" 2)
if (USE_THREAD)

  # N.B. these do NOT have a float type (e.g. DOUBLE) defined!
  GenerateNamedObjects("gemm_thread_m.c;gemm_thread_n.c;gemm_thread_mn.c;gemm_thread_variable.c;syrk_thread.c;level3_thread_lock.c" "" "" 0 "" "" 1)

  if (NOT USE_SIMPLE_THREADED_LEVEL3)
    GenerateCombinationObjects("syrk_k.c" "LOWER;TRANS" "U;N" "THREADED_LEVEL3" 2 "syrk_thread")
    GenerateCombinationObjects("symm_k.c" "RSIDE;LOWER" "L;U" "THREADED_LEVEL3;NN" 2 "symm_thread")
  endif ()
endif ()

foreach (float_type ${FLOAT_TYPES})
  GenerateNamedObjects("gemm_batch_thread.c" "" "gemm_batch_thread" 0 "" "" false ${float_type})

  if (${float_type} STREQUAL "COMPLEX" OR ${float_type} STREQUAL "ZCOMPLEX")
    GenerateCombinationObjects("zherk_kernel.c" "LOWER;CONJ" "U;N" "HERK" 2 "herk_kernel" false ${float_type})
    # TRANS needs to be set/unset when CONJ is set/unset, so can't use it as a combination
    GenerateCombinationObjects("zherk_k.c" "LOWER" "U" "HERK" 3 "herk_N" false ${float_type})
    GenerateCombinationObjects("zherk_k.c" "LOWER" "U" "HERK;TRANS;CONJ" 3 "herk_C" false ${float_type})

    # Need to set CONJ for trmm and trsm
    GenerateCombinationObjects("trmm_L.c" "UPPER;UNIT" "L;N" "CONJ" 0 "trmm_LR" false ${float_type})
    GenerateCombinationObjects("trmm_L.c" "UPPER;UNIT" "L;N" "TRANSA;CONJ" 0 "trmm_LC" false ${float_type})
    GenerateCombinationObjects("trmm_R.c" "UPPER;UNIT" "L;N" "CONJ" 0 "trmm_RR" false ${float_type})
    GenerateCombinationObjects("trmm_R.c" "UPPER;UNIT" "L;N" "TRANSA;CONJ" 0 "trmm_RC" false ${float_type})
    GenerateCombinationObjects("trsm_L.c" "UPPER;UNIT" "L;N" "CONJ" 0 "trsm_LR" false ${float_type})
    GenerateCombinationObjects("trsm_L.c" "UPPER;UNIT" "L;N" "TRANSA;CONJ" 0 "trsm_LC" false ${float_type})
    GenerateCombinationObjects("trsm_R.c" "UPPER;UNIT" "L;N" "CONJ" 0 "trsm_RR" false ${float_type})
    GenerateCombinationObjects("trsm_R.c" "UPPER;UNIT" "L;N" "TRANSA;CONJ" 0 "trsm_RC" false ${float_type})

    #hemm
    GenerateCombinationObjects("zhemm_k.c" "LOWER" "U" "NN" 0 "hemm_L" false ${float_type})
    GenerateCombinationObjects("zhemm_k.c" "LOWER" "U" "NC;RSIDE" 0 "hemm_R" false ${float_type})

    #her2k
    GenerateCombinationObjects("zher2k_kernel.c" "LOWER;CONJ" "U;N" "" 2 "her2k_kernel" false ${float_type})
    GenerateNamedObjects("zher2k_k.c" "HER2K" "her2k_UN" false "" "" false ${float_type})
    GenerateNamedObjects("zher2k_k.c" "HER2K;TRANS;CONJ" "her2k_UC" false "" "" false ${float_type})
    GenerateNamedObjects("zher2k_k.c" "HER2K;LOWER" "her2k_LN" false "" "" false ${float_type})
    GenerateNamedObjects("zher2k_k.c" "HER2K;LOWER;TRANS;CONJ" "her2k_LC" false "" "" false ${float_type})

    if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
      #herk
      GenerateCombinationObjects("zherk_k.c" "LOWER" "U" "HERK;THREADED_LEVEL3" 3 "herk_thread_N" false ${float_type})
      GenerateCombinationObjects("zherk_k.c" "LOWER" "U" "HERK;THREADED_LEVEL3;TRANS;CONJ" 3 "herk_thread_C" false ${float_type})

      #hemm
      GenerateCombinationObjects("zhemm_k.c" "LOWER" "U" "NN;THREADED_LEVEL3" 0 "hemm_thread_L" false ${float_type})
      GenerateCombinationObjects("zhemm_k.c" "LOWER" "U" "NC;RSIDE;THREADED_LEVEL3" 0 "hemm_thread_R" false ${float_type})
      #her2k
      GenerateNamedObjects("zher2k_k.c" "HER2K" "her2k_UN" false "" "" false ${float_type})
      GenerateNamedObjects("zher2k_k.c" "HER2K;TRANS;CONJ" "her2k_UC" false "" "" false ${float_type})
      GenerateNamedObjects("zher2k_k.c" "HER2K;LOWER" "her2k_LN" false "" "" false ${float_type})
      GenerateNamedObjects("zher2k_k.c" "HER2K;LOWER;TRANS;CONJ" "her2k_LC" false "" "" false ${float_type})
    endif()

    # special gemm defines for complex
    foreach (gemm_define ${GEMM_COMPLEX_DEFINES})
      string(TOLOWER ${gemm_define} gemm_define_LC)
      GenerateNamedObjects("gemm.c" "${gemm_define}" "gemm_${gemm_define_LC}" false "" "" false ${float_type})
      if(USE_GEMM3M)
	GenerateNamedObjects("gemm3m.c" "${gemm_define}" "gemm3m_${gemm_define_LC}" false "" "" false ${float_type})
      endif()
      if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
        GenerateNamedObjects("gemm.c" "${gemm_define};THREADED_LEVEL3" "gemm_thread_${gemm_define_LC}" false "" "" false ${float_type})
	if(USE_GEMM3M)
	  GenerateNamedObjects("gemm3m.c" "${gemm_define};THREADED_LEVEL3" "gemm3m_thread_${gemm_define_LC}" false "" "" false ${float_type})
	endif()
      endif ()
    endforeach ()

    if ( BUILD_COMPLEX16 AND NOT  BUILD_DOUBLE)
    foreach (gemm_define ${GEMM_COMPLEX_DEFINES})
      string(TOLOWER ${gemm_define} gemm_define_LC)
      if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
	      GenerateNamedObjects("gemm.c" "${gemm_define};THREADED_LEVEL3" "gemm_thread_${gemm_define_LC}" false "" "" false "DOUBLE" )
    endif()
    endforeach()
    endif ()
    if ( BUILD_COMPLEX AND NOT  BUILD_SINGLE)
    foreach (gemm_define ${GEMM_COMPLEX_DEFINES})
      string(TOLOWER ${gemm_define} gemm_define_LC)
      if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
	      GenerateNamedObjects("gemm.c" "${gemm_define};THREADED_LEVEL3" "gemm_thread_${gemm_define_LC}" false "" "" false "SINGLE" )
    endif()
    endforeach()
    endif ()

    # for gemm3m
    if(USE_GEMM3M)
      foreach (GEMM_DEFINE ${GEMM_DEFINES})
	string(TOLOWER ${GEMM_DEFINE} GEMM_DEFINE_LC)
	GenerateNamedObjects("gemm3m.c" "${GEMM_DEFINE}" "gemm3m_${GEMM_DEFINE_LC}" false "" "" false ${float_type})
	if (USE_THREAD AND NOT USE_SIMPLE_THREADED_LEVEL3)
	  GenerateNamedObjects("gemm3m.c" "${GEMM_DEFINE};THREADED_LEVEL3" "gemm3m_thread_${GEMM_DEFINE_LC}" false "" "" false ${float_type})
	endif ()
      endforeach ()
    endif()
  endif ()
endforeach ()

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

add_library(driver_level3 OBJECT ${OPENBLAS_SRC})

if (USE_OPENMP)
  target_link_libraries(driver_level3 OpenMP::OpenMP_C)
endif()
