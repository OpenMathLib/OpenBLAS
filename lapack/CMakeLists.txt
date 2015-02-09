
include_directories(${CMAKE_SOURCE_DIR})


set(LAPACK_SOURCES
  getf2/getf2_k.c
  getrf/getrf_single.c
  potrf/potrf_U_single.c
  potrf/potrf_L_single.c
  potf2/potf2_U.c
  potf2/potf2_L.c
  lauu2/lauu2_U.c
  lauu2/lauu2_L.c
  lauum/lauum_U_single.c
  lauum/lauum_L_single.c
)

# sources that need TRANS set
set(TRANS_SOURCES
  getrs/getrs_single.c
)

# sources that need UNIT set
set(UNIT_SOURCES
  trtri/trtri_U_single.c
  trtri/trtri_L_single.c
)

set(UNIT_SOURCES2
  trti2/trti2_U.c
  trti2/trti2_L.c
)

# TODO: getrs needs to be compiled with up to TRANS=4 in the complex case
set(ZLAPACK_SOURCES
  getf2/zgetf2_k.c
  getrf/getrf_single.c
  getrs/zgetrs_single.c
  potrf/potrf_U_single.c
  potrf/potrf_L_single.c
  potf2/potf2_U.c
  potf2/potf2_L.c
  lauu2/zlauu2_U.c
  lauu2/zlauu2_L.c
  lauum/lauum_U_single.c
  lauum/lauum_L_single.c
  trti2/ztrti2_U.c
  trti2/ztrti2_L.c
  trtri/trtri_U_single.c
  trtri/trtri_L_single.c
)

GenerateNamedObjects("${LAPACK_SOURCES}" "DOUBLE")

# TODO: laswp needs arch specific code
GenerateNamedObjects("laswp/generic/laswp_k.c" "DOUBLE" "" "laswp_plus")
GenerateNamedObjects("laswp/generic/laswp_k.c" "DOUBLE" "MINUS" "laswp_minus")

if (SMP)

  if (USE_OPENMP)
    set(GETRF_SRC getrf/getrf_parallel_omp.c)
  else ()
    set(GETRF_SRC getrf/getrf_parallel.c)
  endif ()

  set(PARALLEL_SOURCES
    ${GETRF_SRC}
    potrf/potrf_U_parallel.c
    potrf/potrf_L_parallel.c
    lauum/lauum_U_parallel.c
    lauum/lauum_L_parallel.c
  )

  list(APPEND TRANS_SOURCES
    getrs/getrs_parallel.c
  )

  list(APPEND UNIT_SOURCES
    trtri/trtri_U_parallel.c
    trtri/trtri_L_parallel.c
  )

  set(ZPARALLEL_SOURCES
    ${GETRF_SRC}
    getrs/zgetrs_parallel.c
    potrf/potrf_U_parallel.c
    potrf/potrf_L_parallel.c
    lauum/lauum_U_parallel.c
    lauum/lauum_L_parallel.c
    trtri/trtri_U_parallel.c
    trtri/trtri_L_parallel.c
  )

  GenerateNamedObjects("${PARALLEL_SOURCES}" "DOUBLE" "" "" 0)
endif ()

GenerateCombinationObjects("${TRANS_SOURCES}" "TRANS" "N" "DOUBLE" "" 4)
GenerateCombinationObjects("${UNIT_SOURCES}" "UNIT" "N" "DOUBLE" "" 4)
GenerateCombinationObjects("${UNIT_SOURCES2}" "UNIT" "N" "DOUBLE" "" 0)

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

