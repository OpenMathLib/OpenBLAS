
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
GenerateCombinationObjects("${UL_SOURCES}" "LOWER" "U" "" 1 "" "" 3)

# objects that need TRANSA and UNIT set
# N.B. BLAS wants to put the U/L from the filename in the *MIDDLE* because of course why not have a different naming scheme for every single object -hpa
GenerateCombinationObjects("${NU_SOURCES}" "TRANSA;UNIT" "N;N" "" 3 "" "" 3)

# gbmv uses a lowercase n and t. WHY? WHO KNOWS!
GenerateNamedObjects("gbmv_k.c" "" "gbmv_n" false "" "" "" 3)
GenerateNamedObjects("gbmv_k.c" "TRANS" "gbmv_t" false "" "" "" 3)

if (SMP)

  # gbmv uses a lowercase n and t. N.B. this uses TRANSA where gbmv.c uses TRANS. Intentional?
  GenerateNamedObjects("gbmv_thread.c" "" "gbmv_thread_n")
  GenerateNamedObjects("gbmv_thread.c" "TRANSA" "gbmv_thread_t")

  GenerateNamedObjects("gemv_thread.c" "" "gemv_thread_n")
  GenerateNamedObjects("gemv_thread.c" "TRANSA" "gemv_thread_t")

  GenerateNamedObjects("ger_thread.c")

  set(UL_SMP_SOURCES
    symv_thread.c
    syr_thread.c
    syr2_thread.c
    spr_thread.c
    spr2_thread.c
    spmv_thread.c
    sbmv_thread.c
  )

  GenerateCombinationObjects("${UL_SMP_SOURCES}" "LOWER" "U" "" 2)

  set(NU_SMP_SOURCES
    trmv_thread.c
    tpmv_thread.c
    tbmv_thread.c
  )

  GenerateCombinationObjects("${NU_SMP_SOURCES}" "TRANSA;LOWER;UNIT" "N;U;N" "" 2)

endif ()

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

