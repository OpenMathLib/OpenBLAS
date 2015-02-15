
include_directories(${CMAKE_SOURCE_DIR})

# Makeflie

if (${ARCH} STREQUAL "x86")
  GenerateNamedObjects("${KERNELDIR}/cpuid.S" "")
endif ()

# TODO: Read from ${KERNELDIR}/KERNEL - some architectures use a different lsame
set(LSAME_KERNEL lsame.S)
set(SCABS_KERNEL cabs.S)
set(DCABS_KERNEL cabs.S)
# don't use float type name mangling here
GenerateNamedObjects("${KERNELDIR}/${LSAME_KERNEL}" "F_INTERFACE" "lsame" 0 "" "" 1)
GenerateNamedObjects("${KERNELDIR}/${SCABS_KERNEL}" "COMPLEX;F_INTERFACE" "scabs1" "" "" 1)
GenerateNamedObjects("${KERNELDIR}/${DCABS_KERNEL}" "DOUBLE;COMPLEX;F_INTERFACE" "dcabs1" 0 "" "" 1)


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

GenerateNamedObjects("${KERNELDIR}/${DAMAXKERNEL}" "USE_ABS" "amax_k")
GenerateNamedObjects("${KERNELDIR}/${DAMINKERNEL}" "USE_ABS;USE_MIN" "amin_k")
GenerateNamedObjects("${KERNELDIR}/${DMAXKERNEL}" "" "max_k")
GenerateNamedObjects("${KERNELDIR}/${DMINKERNEL}" "" "min_k")
GenerateNamedObjects("${KERNELDIR}/${IDAMAXKERNEL}" "USE_ABS" "i*amax_k")
GenerateNamedObjects("${KERNELDIR}/${IDAMINKERNEL}" "USE_ABS;USE_MIN" "i*amin_k")
GenerateNamedObjects("${KERNELDIR}/${IDMAXKERNEL}" "" "i*max_k")
GenerateNamedObjects("${KERNELDIR}/${IDMINKERNEL}" "" "i*min_k")
GenerateNamedObjects("${KERNELDIR}/${DASUMKERNEL}" "" "asum_k")
GenerateNamedObjects("${KERNELDIR}/${DAXPYKERNEL}" "" "axpy_k")
GenerateNamedObjects("${KERNELDIR}/${DCOPYKERNEL}" "C_INTERFACE" "copy_k")
GenerateNamedObjects("${KERNELDIR}/${DDOTKERNEL}" "" "dot_k")
GenerateNamedObjects("${KERNELDIR}/${DNRM2KERNEL}" "" "nrm2_k")
GenerateNamedObjects("${KERNELDIR}/${DROTKERNEL}" "" "rot_k")
GenerateNamedObjects("${KERNELDIR}/${DSCALKERNEL}" "" "scal_k")
GenerateNamedObjects("${KERNELDIR}/${DSWAPKERNEL}" "" "swap_k")
GenerateNamedObjects("${KERNELDIR}/${DAXPBYKERNEL}" "" "axpby_k")

# Makefile.L2

GenerateNamedObjects("${KERNELDIR}/gemv_n.S" "DOUBLE")
GenerateNamedObjects("${KERNELDIR}/gemv_t.S" "TRANS")
GenerateCombinationObjects("generic/symv_k.c" "LOWER" "U" "" 1)
GenerateNamedObjects("generic/ger.c" "" "ger_k")

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

GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "" "gemm_kernel")

if (DGEMMINCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMINCOPY}" "" "${DGEMMINCOPYOBJ}")
endif ()

if (DGEMMITCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMITCOPY}" "" "${DGEMMITCOPYOBJ}")
endif ()

if (DGEMMONCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMONCOPY}" "" "${DGEMMONCOPYOBJ}")
endif ()

if (DGEMMOTCOPY)
  GenerateNamedObjects("${KERNELDIR}/${DGEMMOTCOPY}" "" "${DGEMMOTCOPYOBJ}")
endif ()

GenerateNamedObjects("${KERNELDIR}/${DGEMM_BETA}")
GenerateCombinationObjects("${KERNELDIR}/${DGEMMKERNEL}" "LEFT;TRANSA" "R;N" "TRMMKERNEL" 2 "trmm_kernel")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "UPPER;LN;TRSMKERNEL" "trsm_kernel_LN")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "LT;TRSMKERNEL" "trsm_kernel_LT")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "UPPER;RN;TRSMKERNEL" "trsm_kernel_RN")
GenerateNamedObjects("${KERNELDIR}/${DGEMMKERNEL}" "RT;TRSMKERNEL" "trsm_kernel_RT")

# These don't use a scheme that is easy to iterate over - the filenames have part of the DEFINE codes in them, for UPPER/TRANS but not for UNIT/OUTER. Also TRANS is not passed in as a define.
# Could simplify it a bit by pairing up by -UUNIT/-DUNIT.
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "UNIT" "trmm_iunucopy")
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "" "trmm_iunncopy")
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "OUTER;UNIT" "trmm_ounucopy")
GenerateNamedObjects("generic/trmm_uncopy_${DGEMM_UNROLL_M}.c" "OUTER" "trmm_ounncopy")

GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "LOWER;UNIT" "trmm_ilnucopy")
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "LOWER" "trmm_ilnncopy")
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trmm_olnucopy")
GenerateNamedObjects("generic/trmm_lncopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER" "trmm_olnncopy")

GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "UNIT" "trmm_iutucopy")
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "" "trmm_iutncopy")
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "OUTER;UNIT" "trmm_outucopy")
GenerateNamedObjects("generic/trmm_utcopy_${DGEMM_UNROLL_M}.c" "OUTER" "trmm_outncopy")

GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "LOWER;UNIT" "trmm_iltucopy")
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "LOWER" "trmm_iltncopy")
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trmm_oltucopy")
GenerateNamedObjects("generic/trmm_ltcopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER" "trmm_oltncopy")

GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "UNIT" "trsm_iunucopy")
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "" "trsm_iunncopy")
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "OUTER;UNIT" "trsm_ounucopy")
GenerateNamedObjects("generic/trsm_uncopy_${DGEMM_UNROLL_M}.c" "OUTER" "trsm_ounncopy")

GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "LOWER;UNIT" "trsm_ilnucopy")
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "LOWER" "trsm_ilnncopy")
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trsm_olnucopy")
GenerateNamedObjects("generic/trsm_lncopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER" "trsm_olnncopy")

GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "UNIT" "trsm_iutucopy")
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "" "trsm_iutncopy")
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "OUTER;UNIT" "trsm_outucopy")
GenerateNamedObjects("generic/trsm_utcopy_${DGEMM_UNROLL_M}.c" "OUTER" "trsm_outncopy")

GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "LOWER;UNIT" "trsm_iltucopy")
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "LOWER" "trsm_iltncopy")
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trsm_oltucopy")
GenerateNamedObjects("generic/trsm_ltcopy_${DGEMM_UNROLL_M}.c" "OUTER;LOWER" "trsm_oltncopy")

GenerateNamedObjects("generic/symm_ucopy_${DGEMM_UNROLL_N}.c" "OUTER" "symm_outcopy")
GenerateNamedObjects("generic/symm_ucopy_${DGEMM_UNROLL_N}.c" "" "symm_iutcopy")

GenerateNamedObjects("generic/symm_lcopy_${DGEMM_UNROLL_N}.c" "LOWER;OUTER" "symm_oltcopy")
GenerateNamedObjects("generic/symm_lcopy_${DGEMM_UNROLL_N}.c" "LOWER" "symm_iltcopy")

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

GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_CN}" "" "domatcopy_k_cn")
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_RN}" "ROWM" "domatcopy_k_rn")
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_CT}" "" "domatcopy_k_ct")
GenerateNamedObjects("${KERNELDIR}/${DOMATCOPY_RT}" "ROWM" "domatcopy_k_rt")

# Makefile.LA
#DBLASOBJS	+= dneg_tcopy$(TSUFFIX).$(SUFFIX) dlaswp_ncopy$(TSUFFIX).$(SUFFIX)

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS
