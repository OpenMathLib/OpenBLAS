###############################################################################
# Copyright (c) 2025, The OpenBLAS Project
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
# 3. Neither the name of the OpenBLAS project nor the names of
#    its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})


set(BLAS1_SOURCES
  copy.c
  nrm2.c
)

set(BLAS1_REAL_ONLY_SOURCES
  rotm.c rotmg.c # N.B. these do not have complex counterparts
  rot.c
  asum.c
  sum.c
)

# these will have 'z' prepended for the complex version
set(BLAS1_MANGLED_SOURCES
  axpy.c swap.c
  scal.c
  dot.c
  rotg.c
  axpby.c
)

# TODO: USE_NETLIB_GEMV should switch gemv.c to netlib/*gemv.f
# these all have 'z' sources for complex versions
set(BLAS2_SOURCES
  gemv.c ger.c
  trsv.c trmv.c 
  syr2.c gbmv.c
  sbmv.c spmv.c
  spr2.c
  tbsv.c tbmv.c
  tpsv.c tpmv.c
)

set(BLAS2_REAL_ONLY_SOURCES
  symv.c syr.c spr.c
)
set(BLAS2_COMPLEX_LAPACK_SOURCES
  symv.c syr.c spr.c
)

set(BLAS2_COMPLEX_ONLY_MANGLED_SOURCES
  hemv.c hbmv.c
  her.c her2.c
  hpmv.c hpr.c
  hpr2.c
)

# these do not have separate 'z' sources
set(BLAS3_SOURCES
  gemm.c symm.c
  trsm.c syrk.c syr2k.c gemmt.c
)

set(BLAS3_MANGLED_SOURCES
  omatcopy.c imatcopy.c
  geadd.c
)

# generate the BLAS objs once with and once without cblas
set (CBLAS_FLAGS "")

if (NOT DEFINED NO_FBLAS)
  list(APPEND CBLAS_FLAGS 0)
endif ()

if (NOT DEFINED NO_CBLAS)
  list(APPEND CBLAS_FLAGS 1)
endif ()

foreach (CBLAS_FLAG ${CBLAS_FLAGS})

  # TODO: don't compile complex sources with cblas for now, the naming schemes are all different and they will have to be handled separately from SINGLE/DOUBLE
  set(DISABLE_COMPLEX 0)
  set(MANGLE_COMPLEX 3)
  if (CBLAS_FLAG EQUAL 1)
#    set(DISABLE_COMPLEX 1)
#    set(MANGLE_COMPLEX 1)
  endif ()
  GenerateNamedObjects("${BLAS1_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false ${DISABLE_COMPLEX})
  GenerateNamedObjects("${BLAS1_REAL_ONLY_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false 1)
  GenerateNamedObjects("${BLAS1_MANGLED_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false ${MANGLE_COMPLEX})
  GenerateNamedObjects("${BLAS2_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false ${MANGLE_COMPLEX})
  GenerateNamedObjects("${BLAS2_REAL_ONLY_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false 1)
  if (NOT DEFINED NO_LAPACK)
  GenerateNamedObjects("${BLAS2_COMPLEX_LAPACK_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false ${MANGLE_COMPLEX})
  endif ()
  GenerateNamedObjects("${BLAS2_COMPLEX_ONLY_MANGLED_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false 4)
  GenerateNamedObjects("${BLAS3_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false ${DISABLE_COMPLEX})
  GenerateNamedObjects("${BLAS3_MANGLED_SOURCES}" "" "" ${CBLAS_FLAG} "" "" false ${MANGLE_COMPLEX})

  GenerateNamedObjects("xerbla.c" "" "xerbla" ${CBLAS_FLAG} "" "" true)
  #sdsdot, dsdot
  if (BUILD_SINGLE OR BUILD_DOUBLE)
  GenerateNamedObjects("sdsdot.c" "" "sdsdot" ${CBLAS_FLAG} "" "" true "SINGLE")
  GenerateNamedObjects("gemm_batch.c" "" "gemm_batch" ${CBLAS_FLAG} "" "" false)
  GenerateNamedObjects("gemm_batch_strided.c" "" "gemm_batch_strided" ${CBLAS_FLAG} "" "" false)
  endif ()
if (BUILD_DOUBLE)
  GenerateNamedObjects("dsdot.c" "" "dsdot" ${CBLAS_FLAG} "" "" true "SINGLE")
endif ()

  # trmm is trsm with a compiler flag set
  GenerateNamedObjects("trsm.c" "TRMM" "trmm" ${CBLAS_FLAG})
  
  # gemmtr is gemmt under the name adopted by the Reference BLAS
  GenerateNamedObjects("gemm.c" "RNAME" "gemmtr" ${CBLAS_FLAG})

  # max and imax are compiled 4 times
  GenerateNamedObjects("max.c" "" "" ${CBLAS_FLAG})
  GenerateNamedObjects("max.c" "USE_ABS" "amax" ${CBLAS_FLAG})
  GenerateNamedObjects("max.c" "USE_ABS;USE_MIN" "amin" ${CBLAS_FLAG})
  GenerateNamedObjects("max.c" "USE_MIN" "min" ${CBLAS_FLAG})

  GenerateNamedObjects("imax.c" "" "i*max" ${CBLAS_FLAG})
  GenerateNamedObjects("imax.c" "USE_ABS" "i*amax" ${CBLAS_FLAG})
  GenerateNamedObjects("imax.c" "USE_ABS;USE_MIN" "i*amin" ${CBLAS_FLAG})
  GenerateNamedObjects("imax.c" "USE_MIN" "i*min" ${CBLAS_FLAG})

if (BUILD_BFLOAT16)
  GenerateNamedObjects("scal.c" "BGEMM" "bscal" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("bf16dot.c" "" "sbdot" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("gemm.c" "BGEMM" "bgemm" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("gemm.c" "" "sbgemm" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("sbgemmt.c" "" "sbgemmt" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("sbgemmt.c" "RNAME" "sbgemmtr" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("sbgemv.c" "BGEMM" "bgemv" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("sbgemv.c" "" "sbgemv" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("tobf16.c" "SINGLE_PREC" "sbstobf16" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("tobf16.c" "DOUBLE_PREC" "sbdtobf16" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("bf16to.c" "SINGLE_PREC" "sbf16tos" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("bf16to.c" "DOUBLE_PREC" "dbf16tod" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("gemm_batch.c" "" "sbgemm_batch" ${CBLAS_FLAG} "" "" true "BFLOAT16")
	GenerateNamedObjects("gemm_batch_strided.c" "" "sbgemm_batch_strided" ${CBLAS_FLAG} "" "" true "BFLOAT16")
endif ()
if (BUILD_HFLOAT16)
  GenerateNamedObjects("gemm.c" "" "shgemm" ${CBLAS_FLAG} "" "" true "HFLOAT16")
endif ()

# complex-specific sources
foreach (float_type ${FLOAT_TYPES})

  if (${float_type} STREQUAL "COMPLEX" OR ${float_type} STREQUAL "ZCOMPLEX")
    GenerateNamedObjects("zaxpy.c" "CONJ" "axpyc" ${CBLAS_FLAG} "" "" false ${float_type})

    GenerateNamedObjects("zger.c" "" "geru" ${CBLAS_FLAG} "" "" false ${float_type})
    GenerateNamedObjects("zger.c" "CONJ" "gerc" ${CBLAS_FLAG} "" "" false ${float_type})
    GenerateNamedObjects("zdot.c" "CONJ" "dotc" ${CBLAS_FLAG} "" "" false ${float_type})
    GenerateNamedObjects("zdot.c" "" "dotu" ${CBLAS_FLAG} "" "" false ${float_type})

    GenerateNamedObjects("symm.c" "HEMM" "hemm" ${CBLAS_FLAG} "" "" false ${float_type})
    GenerateNamedObjects("syrk.c" "HEMM" "herk" ${CBLAS_FLAG} "" "" false ${float_type})
    GenerateNamedObjects("syr2k.c" "HEMM" "her2k" ${CBLAS_FLAG} "" "" false ${float_type})

    if (USE_GEMM3M)
      GenerateNamedObjects("gemm.c" "GEMM3M" "gemm3m" ${CBLAS_FLAG} "" "" false ${float_type})
    endif()
  endif ()
  if (${float_type} STREQUAL "COMPLEX")
    GenerateNamedObjects("zscal.c" "SSCAL" "sscal" ${CBLAS_FLAG} "" "" false "COMPLEX")
    GenerateNamedObjects("nrm2.c" "" "scnrm2" ${CBLAS_FLAG} "" "" true "COMPLEX")
    GenerateNamedObjects("zrot.c" "" "csrot" ${CBLAS_FLAG} "" "" true "COMPLEX")
    GenerateNamedObjects("max.c" "USE_ABS;USE_MIN" "scamin" ${CBLAS_FLAG} "" "" true "COMPLEX")
    GenerateNamedObjects("max.c" "USE_ABS" "scamax" ${CBLAS_FLAG} "" "" true "COMPLEX")
    GenerateNamedObjects("asum.c" "" "scasum" ${CBLAS_FLAG} "" "" true "COMPLEX")
    GenerateNamedObjects("sum.c" "" "scsum" ${CBLAS_FLAG} "" "" true "COMPLEX")
    GenerateNamedObjects("gemm_batch.c" "" "cgemm_batch" ${CBLAS_FLAG} "" "" true "COMPLEX")
    GenerateNamedObjects("gemm_batch_strided.c" "" "cgemm_batch_strided" ${CBLAS_FLAG} "" "" true "COMPLEX")
  endif ()
  if (${float_type} STREQUAL "ZCOMPLEX")
    GenerateNamedObjects("zscal.c" "SSCAL" "dscal" ${CBLAS_FLAG} "" "" false "ZCOMPLEX")
    GenerateNamedObjects("nrm2.c" "" "dznrm2" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
    GenerateNamedObjects("zrot.c" "" "zdrot" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
    GenerateNamedObjects("max.c" "USE_ABS;USE_MIN" "dzamin" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
    GenerateNamedObjects("max.c" "USE_ABS" "dzamax" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
    GenerateNamedObjects("asum.c" "" "dzasum" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
    GenerateNamedObjects("sum.c" "" "dzsum" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
    GenerateNamedObjects("gemm_batch.c" "" "zgemm_batch" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
    GenerateNamedObjects("gemm_batch_strided.c" "" "zgemm_batch_strided" ${CBLAS_FLAG} "" "" true "ZCOMPLEX")
  endif ()
endforeach ()

endforeach ()

#Special functions for CBLAS
if (NOT DEFINED NO_CBLAS)
  foreach (float_type ${FLOAT_TYPES})
  if (${float_type} STREQUAL "COMPLEX" OR ${float_type} STREQUAL "ZCOMPLEX")
    #cblas_dotc_sub cblas_dotu_sub
    GenerateNamedObjects("zdot.c" "FORCE_USE_STACK" "dotu_sub" 1 "" "" false ${float_type})
    GenerateNamedObjects("zdot.c" "FORCE_USE_STACK;CONJ" "dotc_sub" 1 "" "" false ${float_type})
  endif()
  endforeach ()
endif()

if (NOT NO_LAPACK)
  set(LAPACK_SOURCES
    lapack/gesv.c
  )

  # prepend z for complex versions
  set(LAPACK_MANGLED_SOURCES
    lapack/getrf.c lapack/getrs.c lapack/potrf.c lapack/getf2.c
    lapack/potf2.c lapack/laswp.c lapack/lauu2.c
    lapack/lauum.c lapack/trti2.c lapack/trtri.c
  )

  GenerateNamedObjects("${LAPACK_SOURCES}")
  if (NOT RELAPACK_REPLACE)
  GenerateNamedObjects("${LAPACK_MANGLED_SOURCES}" "" "" 0 "" "" 0 3)
  else ()
  GenerateNamedObjects("lapack/getrs.c" "" "" 0 "" "" 0 3)
  GenerateNamedObjects("lapack/getf2.c" "" "" 0 "" "" 0 3)
  GenerateNamedObjects("lapack/potf2.c" "" "" 0 "" "" 0 3)
  GenerateNamedObjects("lapack/laswp.c" "" "" 0 "" "" 0 3)
  GenerateNamedObjects("lapack/lauu2.c" "" "" 0 "" "" 0 3)
  GenerateNamedObjects("lapack/trti2.c" "" "" 0 "" "" 0 3)
  endif()
  GenerateNamedObjects("lapack/laed3.c" "" "" 0 "" "" 0 1)
endif ()

if ( BUILD_COMPLEX AND NOT  BUILD_SINGLE)
	GenerateNamedObjects("scal.c" "" "scal" 0 "" "" false "SINGLE")
	GenerateNamedObjects("copy.c" "" "copy" 0 "" "" false "SINGLE")
	GenerateNamedObjects("dot.c" "" "dot" 0 "" "" false "SINGLE")
	GenerateNamedObjects("rot.c" "" "rot" 0 "" "" false "SINGLE")
	GenerateNamedObjects("nrm2.c" "" "nrm2" 0 "" "" false "SINGLE")
	GenerateNamedObjects("gemv.c" "" "gemv" 0 "" "" false "SINGLE")
	GenerateNamedObjects("gemm.c" "" "gemm" 0 "" "" false "SINGLE")
	GenerateNamedObjects("gemm_batch.c" "" "gemm_batch" 0 "" "" false "SINGLE")
	GenerateNamedObjects("gemm_batch_strided.c" "" "gemm_batch_strided" 0 "" "" false "SINGLE")
	GenerateNamedObjects("asum.c" "" "asum" 0 "" "" false "SINGLE")
	GenerateNamedObjects("swap.c" "" "swap" 0 "" "" false "SINGLE")
	GenerateNamedObjects("axpy.c" "" "axpy" 0 "" "" false "SINGLE")
	GenerateNamedObjects("imax.c" "USE_ABS" "i*amax" 0 "" "" false "SINGLE")
endif ()
if ( BUILD_COMPLEX16 AND NOT  BUILD_DOUBLE)
	GenerateNamedObjects("scal.c" "" "scal" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("copy.c" "" "copy" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("dot.c" "" "dot" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("rot.c" "" "rot" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("nrm2.c" "" "nrm2" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("gemv.c" "" "gemv" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("gemm.c" "" "gemm" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("gemm_batch.c" "" "gemm_batch" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("gemm_batch_strided.c" "" "gemm_batch_strided" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("asum.c" "" "asum" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("swap.c" "" "swap" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("axpy.c" "" "axpy" 0 "" "" false "DOUBLE")
	GenerateNamedObjects("imax.c" "USE_ABS" "i*amax" 0 "" "" false "DOUBLE")
endif ()

if ( BUILD_BFLOAT16 AND NO_FBLAS )
	GenerateNamedObjects("tobf16.c" "SINGLE_PREC" "sbstobf16" 0 "" "" true "BFLOAT16")
	GenerateNamedObjects("bf16to.c" "SINGLE_PREC" "sbf16tos" 0 "" "" true "BFLOAT16")
endif()

add_library(interface OBJECT ${OPENBLAS_SRC})

if (USE_OPENMP)
  target_link_libraries(interface OpenMP::OpenMP_C)
endif()
