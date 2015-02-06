
include_directories(${CMAKE_SOURCE_DIR})

# TODO: need to read ${KERNELDIR}/KERNEL into CMake vars


# Makefile.L1

# these are using hardcoded filenames for now, should get them from the KERNEL vars, e.g. DAMAXKERNEL instead of amax.S
GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "USE_ABS" "amax_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "USE_ABS;USE_MIN" "amin_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "" "max_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/amax.S" "DOUBLE" "" "min_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "USE_ABS" "i*amax_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "USE_ABS;USE_MIN" "i*amin_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "" "i*max_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

GenerateNamedObjects("${KERNELDIR}/iamax.S" "DOUBLE" "" "i*min_k" 0)
list(APPEND DBLAS_OBJS ${OBJ_LIST_OUT})

#DBLASOBJS	+= \
#	dasum_k$(TSUFFIX).$(SUFFIX) daxpy_k$(TSUFFIX).$(SUFFIX) dcopy_k$(TSUFFIX).$(SUFFIX) ddot_k$(TSUFFIX).$(SUFFIX) \
#	dnrm2_k$(TSUFFIX).$(SUFFIX) drot_k$(TSUFFIX).$(SUFFIX) dscal_k$(TSUFFIX).$(SUFFIX) dswap_k$(TSUFFIX).$(SUFFIX) \
#	daxpby_k$(TSUFFIX).$(SUFFIX)

# Makefile.L2
#DBLASOBJS	+= \
#	dgemv_n$(TSUFFIX).$(SUFFIX) dgemv_t$(TSUFFIX).$(SUFFIX) dsymv_U$(TSUFFIX).$(SUFFIX) dsymv_L$(TSUFFIX).$(SUFFIX) \
#	dger_k$(TSUFFIX).$(SUFFIX)


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
