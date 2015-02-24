
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

# sources that need to be compiled several times, for UNIT, TRANSA
set(U_SOURCES
  trmv_U.c
  tbmv_U.c
  tbsv_U.c
  tpmv_U.c
  tpsv_U.c
  trsv_U.c
)

set(L_SOURCES
  trmv_L.c
  tbmv_L.c
  tbsv_L.c
  tpmv_L.c
  tpsv_L.c
  trsv_L.c
)

set(UL_SMP_SOURCES
  symv_thread.c
  syr_thread.c
  syr2_thread.c
  spr_thread.c
  spr2_thread.c
  spmv_thread.c
  sbmv_thread.c
)

set(NU_SMP_SOURCES
  trmv_thread.c
  tpmv_thread.c
  tbmv_thread.c
)

# objects that need LOWER set
GenerateCombinationObjects("${UL_SOURCES}" "LOWER" "U" "" 1 "" "" 3)

# gbmv uses a lowercase n and t
GenerateNamedObjects("gbmv_k.c" "" "gbmv_n" false "" "" "" 3)
GenerateNamedObjects("gbmv_k.c" "TRANS" "gbmv_t" false "" "" "" 3)

# special defines for complex
foreach (float_type ${FLOAT_TYPES})

  if (SMP)
    GenerateNamedObjects("gemv_thread.c" "" "gemv_thread_n" false "" "" false ${float_type})
    GenerateNamedObjects("gemv_thread.c" "TRANSA" "gemv_thread_t" false "" "" false ${float_type})

    GenerateNamedObjects("gbmv_thread.c" "" "gbmv_thread_n" false "" "" false ${float_type})
    GenerateNamedObjects("gbmv_thread.c" "TRANSA" "gbmv_thread_t" false "" "" false ${float_type})
  endif ()

  if (${float_type} STREQUAL "COMPLEX" OR ${float_type} STREQUAL "ZCOMPLEX")

    foreach (u_source ${U_SOURCES})
      string(REGEX MATCH "[a-z]+" op_name ${u_source})
      GenerateCombinationObjects("z${u_source}" "UNIT" "N" "TRANSA=1" 0 "${op_name}_NU" false ${float_type})
      GenerateCombinationObjects("z${u_source}" "UNIT" "N" "TRANSA=2" 0 "${op_name}_TL" false ${float_type})
      GenerateCombinationObjects("z${u_source}" "UNIT" "N" "TRANSA=3" 0 "${op_name}_RU" false ${float_type})
      GenerateCombinationObjects("z${u_source}" "UNIT" "N" "TRANSA=4" 0 "${op_name}_CL" false ${float_type})
    endforeach ()

    foreach (l_source ${L_SOURCES})
      string(REGEX MATCH "[a-z]+" op_name ${l_source})
      GenerateCombinationObjects("z${l_source}" "UNIT" "N" "TRANSA=1" 0 "${op_name}_NL" false ${float_type})
      GenerateCombinationObjects("z${l_source}" "UNIT" "N" "TRANSA=2" 0 "${op_name}_TU" false ${float_type})
      GenerateCombinationObjects("z${l_source}" "UNIT" "N" "TRANSA=3" 0 "${op_name}_RL" false ${float_type})
      GenerateCombinationObjects("z${l_source}" "UNIT" "N" "TRANSA=4" 0 "${op_name}_CU" false ${float_type})
    endforeach ()

    if (SMP)

      GenerateCombinationObjects("trmv_thread.c" "LOWER;UNIT" "U;N" "TRANSA=3" 0 "trmv_thread_R" false ${float_type})
      GenerateCombinationObjects("trmv_thread.c" "LOWER;UNIT" "U;N" "TRANSA=4" 0 "trmv_thread_C" false ${float_type})

      GenerateNamedObjects("gemv_thread.c" "CONJ" "gemv_thread_r" false "" "" false ${float_type})
      GenerateNamedObjects("gemv_thread.c" "CONJ;TRANSA" "gemv_thread_c" false "" "" false ${float_type})
      GenerateNamedObjects("gemv_thread.c" "XCONJ" "gemv_thread_o" false "" "" false ${float_type})
      GenerateNamedObjects("gemv_thread.c" "XCONJ;TRANSA" "gemv_thread_u" false "" "" false ${float_type})
      GenerateNamedObjects("gemv_thread.c" "XCONJ;CONJ" "gemv_thread_s" false "" "" false ${float_type})
      GenerateNamedObjects("gemv_thread.c" "XCONJ;CONJ;TRANSA" "gemv_thread_d" false "" "" false ${float_type})

      GenerateNamedObjects("gbmv_thread.c" "CONJ" "gbmv_thread_r" false "" "" false ${float_type})
      GenerateNamedObjects("gbmv_thread.c" "CONJ;TRANSA" "gbmv_thread_c" false "" "" false ${float_type})
      GenerateNamedObjects("gbmv_thread.c" "XCONJ" "gbmv_thread_o" false "" "" false ${float_type})
      GenerateNamedObjects("gbmv_thread.c" "XCONJ;TRANSA" "gbmv_thread_u" false "" "" false ${float_type})
      GenerateNamedObjects("gbmv_thread.c" "XCONJ;CONJ" "gbmv_thread_s" false "" "" false ${float_type})
      GenerateNamedObjects("gbmv_thread.c" "XCONJ;CONJ;TRANSA" "gbmv_thread_d" false "" "" false ${float_type})

      GenerateNamedObjects("ger_thread.c" "" "ger_thread_U" false "" "" false ${float_type})
      GenerateNamedObjects("ger_thread.c" "CONJ" "ger_thread_C" false "" "" false ${float_type})
      GenerateNamedObjects("ger_thread.c" "XCONJ" "ger_thread_V" false "" "" false ${float_type})
      GenerateNamedObjects("ger_thread.c" "XCONJ;CONJ" "ger_thread_D" false "" "" false ${float_type})

      foreach (nu_smp_src ${NU_SMP_SOURCSE})
        string(REGEX MATCH "[a-z]+" op_name ${nu_smp_src})
        GenerateCombinationObjects("${nu_smp_src}" "LOWER;UNIT" "U;N" "TRANSA=1" 0 "${op_name}_N" false ${float_type})
        GenerateCombinationObjects("${nu_smp_src}" "LOWER;UNIT" "U;N" "TRANSA=2" 0 "${op_name}_T" false ${float_type})
        GenerateCombinationObjects("${nu_smp_src}" "LOWER;UNIT" "U;N" "TRANSA=3" 0 "${op_name}_R" false ${float_type})
        GenerateCombinationObjects("${nu_smp_src}" "LOWER;UNIT" "U;N" "TRANSA=4" 0 "${op_name}_C" false ${float_type})
      endforeach ()
    endif ()

  else ()
    # N.B. BLAS wants to put the U/L from the filename in the *MIDDLE*
    GenerateCombinationObjects("${U_SOURCES};${L_SOURCES}" "TRANSA;UNIT" "N;N" "" 3 "" false ${float_type})
    if (SMP)
      GenerateNamedObjects("ger_thread.c" "" "" false "" "" false ${float_type})
      GenerateCombinationObjects("${NU_SMP_SOURCES}" "TRANSA;LOWER;UNIT" "N;U;N" "" 2 "" false ${float_type})
    endif ()
  endif ()
endforeach ()

if (SMP)
  GenerateCombinationObjects("${UL_SMP_SOURCES}" "LOWER" "U" "" 2)
endif ()

add_library(driver_level2 OBJECT ${OPENBLAS_SRC})
