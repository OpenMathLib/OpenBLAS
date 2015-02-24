
include_directories(${CMAKE_SOURCE_DIR})


set(LAPACK_SOURCES
  getrf/getrf_single.c
  potrf/potrf_U_single.c
  potrf/potrf_L_single.c
  lauum/lauum_U_single.c
  lauum/lauum_L_single.c
)

# add a 'z' to filename for complex version
set(LAPACK_MANGLED_SOURCES
  getf2/getf2_k.c
  lauu2/lauu2_U.c
  lauu2/lauu2_L.c
  potf2/potf2_U.c
  potf2/potf2_L.c
)

# sources that need TRANS set
# this has a 'z' version
set(TRANS_SOURCES
  getrs/getrs_single.c
)

# sources that need UNIT set
# these do NOT have a z version
set(UNIT_SOURCES
  trtri/trtri_U_single.c
  trtri/trtri_L_single.c
)

# these have a 'z' version
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

GenerateNamedObjects("${LAPACK_SOURCES}")
GenerateNamedObjects("${LAPACK_MANGLED_SOURCES}" "" "" false "" "" false 3)

# TODO: laswp needs arch specific code
GenerateNamedObjects("laswp/generic/laswp_k.c" "" "laswp_plus" false "" "" false 3)
GenerateNamedObjects("laswp/generic/laswp_k.c" "MINUS" "laswp_minus" false "" "" false 3)

if (SMP)

  if (USE_OPENMP)
    set(GETRF_SRC getrf/getrf_parallel_omp.c)
  else ()
    set(GETRF_SRC getrf/getrf_parallel.c)
  endif ()

  # these do not have 'z' versions
  set(PARALLEL_SOURCES
    ${GETRF_SRC}
    lauum/lauum_U_parallel.c
    lauum/lauum_L_parallel.c
    potrf/potrf_U_parallel.c
    potrf/potrf_L_parallel.c
  )

  # this has a z version
  list(APPEND TRANS_SOURCES
    getrs/getrs_parallel.c
  )

  # these do NOT have a z version
  list(APPEND UNIT_SOURCES
    trtri/trtri_U_parallel.c
    trtri/trtri_L_parallel.c
  )

  GenerateNamedObjects("${PARALLEL_SOURCES}")
endif ()

GenerateCombinationObjects("${TRANS_SOURCES}" "TRANS" "N" "" 4 "" "" 3)
GenerateCombinationObjects("${UNIT_SOURCES}" "UNIT" "N" "" 4)
GenerateCombinationObjects("${UNIT_SOURCES2}" "UNIT" "N" "" 0 "" "" 3)

add_library(lapack OBJECT ${OPENBLAS_SRC})

