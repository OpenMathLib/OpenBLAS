
include_directories(${CMAKE_SOURCE_DIR})

# TODO: need to read ${KERNELDIR}/KERNEL into CMake vars


# Makefile.L1

# these are using hardcoded filenames for now, should get them from the KERNEL vars, e.g. DAMAXKERNEL instead of amax.S
GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "USE_ABS" "amax_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "USE_ABS;USE_MIN" "amin_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "" "max_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "" "min_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "USE_ABS" "i*amax_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "USE_ABS;USE_MIN" "i*amin_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "" "i*max_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "" "i*min_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/asum.S" "DOUBLE" "" "asum_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/axpy.S" "DOUBLE" "" "axpy_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/copy.S" "DOUBLE" "C_INTERFACE" "copy_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/dot.S" "DOUBLE" "" "dot_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/nrm2.S" "DOUBLE" "" "nrm2_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/rot.S" "DOUBLE" "" "rot_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/scal.S" "DOUBLE" "" "scal_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/swap.S" "DOUBLE" "" "swap_k")
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("arm/axpby.c" "DOUBLE" "" "axpby_k")
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
#DKERNELOBJS	+= \
#	dgemm_kernel$(TSUFFIX).$(SUFFIX) \
#	$(DGEMMINCOPYOBJ) $(DGEMMITCOPYOBJ) \
#	$(DGEMMONCOPYOBJ) $(DGEMMOTCOPYOBJ)

#DBLASOBJS	+= \
#	dgemm_beta$(TSUFFIX).$(SUFFIX) \
#	dtrmm_kernel_LN$(TSUFFIX).$(SUFFIX) dtrmm_kernel_LT$(TSUFFIX).$(SUFFIX) \
#	dtrmm_kernel_RN$(TSUFFIX).$(SUFFIX) dtrmm_kernel_RT$(TSUFFIX).$(SUFFIX) \
#	dtrsm_kernel_LN$(TSUFFIX).$(SUFFIX) dtrsm_kernel_LT$(TSUFFIX).$(SUFFIX) \
#	dtrsm_kernel_RN$(TSUFFIX).$(SUFFIX) dtrsm_kernel_RT$(TSUFFIX).$(SUFFIX) \

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
