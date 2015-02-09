
include_directories(${CMAKE_SOURCE_DIR})

# Makeflie

if (${ARCH} STREQUAL "x86")
  GenerateNamedObjects("${KERNELDIR}/cpuid.S" "")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
endif ()


# Makefile.L1

# TODO: need to read ${KERNELDIR}/KERNEL into CMake vars
set(DAMAXKERNEL amax.S)
set(DAMINKERNEL amax.S)
set(DMAXKERNEL amax.S)
set(DMINKERNEL amax.S)
set(IDAMAXKERNEL iamax.S)
set(IDAMINKERNEL iamax.S)
set(IDMAXKERNEL iamax.S)
set(IDMINKERNEL iamax.S)
set(DASUMKERNEL asum.S)
set(DAXPYKERNEL axpy.S)
set(DCOPYKERNEL copy.S)
set(DDOTKERNEL dot.S)
set(DNRM2KERNEL nrm2.S)
set(DROTKERNEL rot.S)
set(DSCALKERNEL scal.S)
set(DSWAPKERNEL swap.S)
set(DAXPBYKERNEL ../arm/axpby.c)

GenerateNamedObjects("${KERNELDIR}/${DAMAXKERNEL}" "DOUBLE" "USE_ABS" "amax_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DAMINKERNEL}" "DOUBLE" "USE_ABS;USE_MIN" "amin_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DMAXKERNEL}" "DOUBLE" "" "max_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DMINKERNEL}" "DOUBLE" "" "min_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${IDAMAXKERNEL}" "DOUBLE" "USE_ABS" "i*amax_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${IDAMINKERNEL}" "DOUBLE" "USE_ABS;USE_MIN" "i*amin_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${IDMAXKERNEL}" "DOUBLE" "" "i*max_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${IDMINKERNEL}" "DOUBLE" "" "i*min_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DASUMKERNEL}" "DOUBLE" "" "asum_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DAXPYKERNEL}" "DOUBLE" "" "axpy_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DCOPYKERNEL}" "DOUBLE" "C_INTERFACE" "copy_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DDOTKERNEL}" "DOUBLE" "" "dot_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DNRM2KERNEL}" "DOUBLE" "" "nrm2_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DROTKERNEL}" "DOUBLE" "" "rot_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DSCALKERNEL}" "DOUBLE" "" "scal_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DSWAPKERNEL}" "DOUBLE" "" "swap_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/${DAXPBYKERNEL}" "DOUBLE" "" "axpby_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

# Makefile.L2

GenerateNamedObjects("${KERNELDIR}/gemv_n.S" "DOUBLE")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/gemv_t.S" "DOUBLE" "TRANS")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateCombinationObjects("generic/symv_k.c" "LOWER" "U" "DOUBLE" "" 1)
list(APPEND DBLAS_OBJS ${COMBO_OBJ_LIST_OUT})

GenerateNamedObjects("generic/ger.c" "DOUBLE" "" "ger_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

# Makefile.L3

set(DGEMM_BETA ../generic/gemm_beta.c)
set(DGEMMKERNEL gemm_kernel_2x4_penryn.S)

GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

if (DGEMMINCOPY)
  GenerateNamedObjects(${DGEMMINCOPY} "DOUBLE")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
endif ()

if (DGEMMITCOPY)
  GenerateNamedObjects(${DGEMMITCOPY} "DOUBLE")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
endif ()

if (DGEMMONCOPY)
  GenerateNamedObjects(${DGEMMONCOPY} "DOUBLE")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
endif ()

if (DGEMMOTCOPY)
  GenerateNamedObjects(${DGEMMOTCOPY} "DOUBLE")
  list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
endif ()

GenerateNamedObjects("${KERNELDIR}/${DGEMM_BETA}" "DOUBLE")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateCombinationObjects("${KERNELDIR}/${DGEMMKERNEL}" "LEFT;TRANSA" "R;N" "DOUBLE" "TRMMKERNEL" 2 "trmm_kernel")
list(APPEND DBLAS_OBJS ${COMBO_OBJ_LIST_OUT})
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "UPPER;LN;TRSMKERNEL" "trsm_kernel_LN")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "LT;TRSMKERNEL" "trsm_kernel_LT")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "UPPER;RN;TRSMKERNEL" "trsm_kernel_RN")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "RT;TRSMKERNEL" "trsm_kernel_RT")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

# These don't use a scheme that is easy to iterate over - the filenames have part of the DEFINE codes in them, for UPPER/TRANS but not for UNIT/OUTER. Also TRANS is not passed in as a define.
# Could simplify it a bit by pairing up by -UUNIT/-DUNIT.
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trmm_iunucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trmm_iunncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trmm_ounucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trmm_ounncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trmm_ilnucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trmm_ilnncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trmm_olnucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trmm_olnncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trmm_iutucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trmm_iutncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trmm_outucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trmm_outncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trmm_iltucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trmm_iltncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trmm_oltucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trmm_oltncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trsm_iunucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trsm_iunncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trsm_ounucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trsm_ounncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trsm_ilnucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trsm_ilnncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trsm_olnucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trsm_olnncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trsm_iutucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trsm_iutncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trsm_outucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trsm_outncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trsm_iltucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trsm_iltncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trsm_oltucopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trsm_oltncopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/symm_ucopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "OUTER" "symm_outcopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/symm_ucopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "" "symm_iutcopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("generic/symm_lcopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "LOWER;OUTER" "symm_oltcopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("generic/symm_lcopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "LOWER" "symm_iltcopy")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

if (NOT DEFINED DOMATCOPY_CN)
  set(DOMATCOPY_CN ../arm/omatcopy_cn.c)
endif ()
if (NOT DEFINED DOMATCOPY_RN)
  set(DOMATCOPY_RN ../arm/omatcopy_rn.c)
endif ()
if (NOT DEFINED DOMATCOPY_CT)
  set(DOMATCOPY_CT ../arm/omatcopy_ct.c)
endif ()
if (NOT DEFINED DOMATCOPY_RT)
  set(DOMATCOPY_RT ../arm/omatcopy_rt.c)
endif ()

GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_CN}" "DOUBLE" "" "domatcopy_k_cn")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_RN}" "DOUBLE" "ROWM" "domatcopy_k_rn")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_CT}" "DOUBLE" "" "domatcopy_k_ct")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_RT}" "DOUBLE" "ROWM" "domatcopy_k_rt")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

# Makefile.LA
#DBLASOBJS	+= dneg_tcopy$(TSUFFIX).$(SUFFIX) dlaswp_ncopy$(TSUFFIX).$(SUFFIX)

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS
