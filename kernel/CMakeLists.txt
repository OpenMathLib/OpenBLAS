
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

#DBLASOBJS += \
#	dtrmm_iunucopy$(TSUFFIX).$(SUFFIX) dtrmm_iunncopy$(TSUFFIX).$(SUFFIX) \
#	dtrmm_ilnucopy$(TSUFFIX).$(SUFFIX) dtrmm_ilnncopy$(TSUFFIX).$(SUFFIX) \
#	dtrmm_iutucopy$(TSUFFIX).$(SUFFIX) dtrmm_iutncopy$(TSUFFIX).$(SUFFIX) \
#	dtrmm_iltucopy$(TSUFFIX).$(SUFFIX) dtrmm_iltncopy$(TSUFFIX).$(SUFFIX) \
#	dtrmm_ounucopy$(TSUFFIX).$(SUFFIX) dtrmm_ounncopy$(TSUFFIX).$(SUFFIX) \
#	dtrmm_olnucopy$(TSUFFIX).$(SUFFIX) dtrmm_olnncopy$(TSUFFIX).$(SUFFIX) \
#	dtrmm_outucopy$(TSUFFIX).$(SUFFIX) dtrmm_outncopy$(TSUFFIX).$(SUFFIX) \
#	dtrmm_oltucopy$(TSUFFIX).$(SUFFIX) dtrmm_oltncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_iunucopy$(TSUFFIX).$(SUFFIX) dtrsm_iunncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_ilnucopy$(TSUFFIX).$(SUFFIX) dtrsm_ilnncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_iutucopy$(TSUFFIX).$(SUFFIX) dtrsm_iutncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_iltucopy$(TSUFFIX).$(SUFFIX) dtrsm_iltncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_ounucopy$(TSUFFIX).$(SUFFIX) dtrsm_ounncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_olnucopy$(TSUFFIX).$(SUFFIX) dtrsm_olnncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_outucopy$(TSUFFIX).$(SUFFIX) dtrsm_outncopy$(TSUFFIX).$(SUFFIX) \
#	dtrsm_oltucopy$(TSUFFIX).$(SUFFIX) dtrsm_oltncopy$(TSUFFIX).$(SUFFIX) \
#	dsymm_iutcopy$(TSUFFIX).$(SUFFIX) dsymm_iltcopy$(TSUFFIX).$(SUFFIX) \
#	dsymm_outcopy$(TSUFFIX).$(SUFFIX) dsymm_oltcopy$(TSUFFIX).$(SUFFIX)

#DBLASOBJS += \
#	domatcopy_k_cn$(TSUFFIX).$(SUFFIX) domatcopy_k_rn$(TSUFFIX).$(SUFFIX) \
#	domatcopy_k_ct$(TSUFFIX).$(SUFFIX) domatcopy_k_rt$(TSUFFIX).$(SUFFIX)

# Makefile.LA
#DBLASOBJS	+= dneg_tcopy$(TSUFFIX).$(SUFFIX) dlaswp_ncopy$(TSUFFIX).$(SUFFIX)

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS
