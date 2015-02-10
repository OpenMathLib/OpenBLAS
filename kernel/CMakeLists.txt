
include_directories(${CMAKE_SOURCE_DIR})

# Makeflie

if (${ARCH} STREQUAL "x86")
  GenerateNamedObjects("${KERNELDIR}/cpuid.S" "")
endif ()

# TODO: Read from ${KERNELDIR}/KERNEL - some architectures use a different lsame
set(LSAME_KERNEL lsame.S)
set(SCABS_KERNEL cabs.S)
set(DCABS_KERNEL cabs.S)
GenerateNamedObjects("${KERNELDIR}/${LSAME_KERNEL}" "" "F_INTERFACE" "lsame")
GenerateNamedObjects("${KERNELDIR}/${SCABS_KERNEL}" "SINGLE" "COMPLEX;F_INTERFACE" "cabs1")
GenerateNamedObjects("${KERNELDIR}/${DCABS_KERNEL}" "DOUBLE" "COMPLEX;F_INTERFACE" "cabs1")


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
GenerateNamedObjects("${KERNELDIR}/${DAMINKERNEL}" "DOUBLE" "USE_ABS;USE_MIN" "amin_k")
GenerateNamedObjects("${KERNELDIR}/${DMAXKERNEL}" "DOUBLE" "" "max_k")
GenerateNamedObjects("${KERNELDIR}/${DMINKERNEL}" "DOUBLE" "" "min_k")
GenerateNamedObjects("${KERNELDIR}/${IDAMAXKERNEL}" "DOUBLE" "USE_ABS" "i*amax_k")
GenerateNamedObjects("${KERNELDIR}/${IDAMINKERNEL}" "DOUBLE" "USE_ABS;USE_MIN" "i*amin_k")
GenerateNamedObjects("${KERNELDIR}/${IDMAXKERNEL}" "DOUBLE" "" "i*max_k")
GenerateNamedObjects("${KERNELDIR}/${IDMINKERNEL}" "DOUBLE" "" "i*min_k")
GenerateNamedObjects("${KERNELDIR}/${DASUMKERNEL}" "DOUBLE" "" "asum_k")
GenerateNamedObjects("${KERNELDIR}/${DAXPYKERNEL}" "DOUBLE" "" "axpy_k")
GenerateNamedObjects("${KERNELDIR}/${DCOPYKERNEL}" "DOUBLE" "C_INTERFACE" "copy_k")
GenerateNamedObjects("${KERNELDIR}/${DDOTKERNEL}" "DOUBLE" "" "dot_k")
GenerateNamedObjects("${KERNELDIR}/${DNRM2KERNEL}" "DOUBLE" "" "nrm2_k")
GenerateNamedObjects("${KERNELDIR}/${DROTKERNEL}" "DOUBLE" "" "rot_k")
GenerateNamedObjects("${KERNELDIR}/${DSCALKERNEL}" "DOUBLE" "" "scal_k")
GenerateNamedObjects("${KERNELDIR}/${DSWAPKERNEL}" "DOUBLE" "" "swap_k")
GenerateNamedObjects("${KERNELDIR}/${DAXPBYKERNEL}" "DOUBLE" "" "axpby_k")

# Makefile.L2

GenerateNamedObjects("${KERNELDIR}/gemv_n.S" "DOUBLE")
GenerateNamedObjects("${KERNELDIR}/gemv_t.S" "DOUBLE" "TRANS")
GenerateCombinationObjects("generic/symv_k.c" "LOWER" "U" "DOUBLE" "" 1)
GenerateNamedObjects("generic/ger.c" "DOUBLE" "" "ger_k")

# Makefile.L3

# TODO: these are from KERNEL.PENRYN - they should be read in from the appropriate ${KERNELDIR}/KERNEL file
set(DGEMM_BETA ../generic/gemm_beta.c)
set(DGEMMKERNEL gemm_kernel_2x4_penryn.S)
set(DGEMMINCOPY gemm_ncopy_2.S)
set(DGEMMITCOPY gemm_tcopy_2.S)
set(DGEMMONCOPY ../generic/gemm_ncopy_4.c)
set(DGEMMOTCOPY ../generic/gemm_tcopy_4.c)
set(DGEMMINCOPYOBJ gemm_incopy)
set(DGEMMITCOPYOBJ gemm_itcopy)
set(DGEMMONCOPYOBJ gemm_oncopy)
set(DGEMMOTCOPYOBJ gemm_otcopy)

GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "" "gemm_kernel")

if (DGEMMINCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMINCOPY}" "DOUBLE" "" "${DGEMMINCOPYOBJ}")
endif ()

if (DGEMMITCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMITCOPY}" "DOUBLE" "" "${DGEMMITCOPYOBJ}")
endif ()

if (DGEMMONCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMONCOPY}" "DOUBLE" "" "${DGEMMONCOPYOBJ}")
endif ()

if (DGEMMOTCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMOTCOPY}" "DOUBLE" "" "${DGEMMOTCOPYOBJ}")
endif ()

GenerateNamedObjects("${KERNELDIR}/${DGEMM_BETA}" "DOUBLE")
GenerateCombinationObjects("${KERNELDIR}/${DGEMMKERNEL}" "LEFT;TRANSA" "R;N" "DOUBLE" "TRMMKERNEL" 2 "trmm_kernel")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "UPPER;LN;TRSMKERNEL" "trsm_kernel_LN")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "LT;TRSMKERNEL" "trsm_kernel_LT")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "UPPER;RN;TRSMKERNEL" "trsm_kernel_RN")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "DOUBLE" "RT;TRSMKERNEL" "trsm_kernel_RT")

# These don't use a scheme that is easy to iterate over - the filenames have part of the DEFINE codes in them, for UPPER/TRANS but not for UNIT/OUTER. Also TRANS is not passed in as a define.
# Could simplify it a bit by pairing up by -UUNIT/-DUNIT.
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trmm_iunucopy")
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trmm_iunncopy")
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trmm_ounucopy")
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trmm_ounncopy")

GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trmm_ilnucopy")
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trmm_ilnncopy")
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trmm_olnucopy")
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trmm_olnncopy")

GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trmm_iutucopy")
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trmm_iutncopy")
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trmm_outucopy")
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trmm_outncopy")

GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trmm_iltucopy")
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trmm_iltncopy")
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trmm_oltucopy")
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trmm_oltncopy")

GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trsm_iunucopy")
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trsm_iunncopy")
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trsm_ounucopy")
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trsm_ounncopy")

GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trsm_ilnucopy")
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trsm_ilnncopy")
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trsm_olnucopy")
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trsm_olnncopy")

GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "UNIT" "trsm_iutucopy")
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "" "trsm_iutncopy")
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;UNIT" "trsm_outucopy")
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER" "trsm_outncopy")

GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER;UNIT" "trsm_iltucopy")
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "LOWER" "trsm_iltncopy")
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER;UNIT" "trsm_oltucopy")
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "DOUBLE" "OUTER;LOWER" "trsm_oltncopy")

GenerateNamedObjects("generic/symm_ucopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "OUTER" "symm_outcopy")
GenerateNamedObjects("generic/symm_ucopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "" "symm_iutcopy")

GenerateNamedObjects("generic/symm_lcopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "LOWER;OUTER" "symm_oltcopy")
GenerateNamedObjects("generic/symm_lcopy_${DGEMM_UNROLL_N}.c" "DOUBLE" "LOWER" "symm_iltcopy")

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
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_RN}" "DOUBLE" "ROWM" "domatcopy_k_rn")
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_CT}" "DOUBLE" "" "domatcopy_k_ct")
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_RT}" "DOUBLE" "ROWM" "domatcopy_k_rt")

# Makefile.LA
#DBLASOBJS	+= dneg_tcopy$(TSUFFIX).$(SUFFIX) dlaswp_ncopy$(TSUFFIX).$(SUFFIX)

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS
