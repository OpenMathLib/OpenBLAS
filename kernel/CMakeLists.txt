
include_directories(${CMAKE_SOURCE_DIR})
include("${CMAKE_SOURCE_DIR}/cmake/kernel.cmake")

# Makefile

if (DEFINED TARGET_CORE)
  #override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE)
  set(BUILD_KERNEL 1)
  set(KDIR "")
  set(TSUFFIX "_${TARGET_CORE}")
else ()
  set(TARGET_CORE ${CORE})
  set(KDIR "")
  set(TSUFFIX "")
endif ()

SetDefaultL1()
#-include $(KERNELDIR)/KERNEL.$(TARGET_CORE)
#include $(KERNELDIR)/KERNEL
ParseMakefileVars("${KERNELDIR}/KERNEL.${TARGET_CORE}")
ParseMakefileVars("${KERNELDIR}/KERNEL")

if (${ARCH} STREQUAL "x86")
  GenerateNamedObjects("${KERNELDIR}/cpuid.S" "")
endif ()

# don't use float type name mangling here
GenerateNamedObjects("${KERNELDIR}/${LSAME_KERNEL}" "F_INTERFACE" "lsame" false "" "" true)
GenerateNamedObjects("${KERNELDIR}/${SCABS_KERNEL}" "COMPLEX;F_INTERFACE" "scabs1" false "" "" true)
GenerateNamedObjects("${KERNELDIR}/${DCABS_KERNEL}" "DOUBLE;COMPLEX;F_INTERFACE" "dcabs1" false "" "" true)

# Makefile.L1
foreach (float_type ${FLOAT_TYPES})
  # a bit of metaprogramming here to pull out the appropriate KERNEL var
  string(SUBSTRING ${float_type} 0 1 float_char)
  GenerateNamedObjects("${KERNELDIR}/${${float_char}AMAXKERNEL}" "USE_ABS" "amax_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}AMINKERNEL}" "USE_ABS;USE_MIN" "amin_k" false "" "" false ${float_type})
  if (DEFINED ${float_char}MAXKERNEL)
    GenerateNamedObjects("${KERNELDIR}/${${float_char}MAXKERNEL}" "" "max_k" false "" "" false ${float_type})
  endif ()
  if (DEFINED ${float_char}MINKERNEL)
    GenerateNamedObjects("${KERNELDIR}/${${float_char}MINKERNEL}" "" "min_k" false "" "" false ${float_type})
  endif ()
  GenerateNamedObjects("${KERNELDIR}/${I${float_char}AMAXKERNEL}" "USE_ABS" "i*amax_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${I${float_char}AMINKERNEL}" "USE_ABS;USE_MIN" "i*amin_k" false "" "" false ${float_type})
  if (DEFINED I${float_char}MAXKERNEL)
    GenerateNamedObjects("${KERNELDIR}/${I${float_char}MAXKERNEL}" "" "i*max_k" false "" "" false ${float_type})
  endif ()
  if (DEFINED I${float_char}MINKERNEL)
    GenerateNamedObjects("${KERNELDIR}/${I${float_char}MINKERNEL}" "" "i*min_k" false "" "" false ${float_type})
  endif ()
  GenerateNamedObjects("${KERNELDIR}/${${float_char}ASUMKERNEL}" "" "asum_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}AXPYKERNEL}" "" "axpy_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}COPYKERNEL}" "C_INTERFACE" "copy_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}DOTKERNEL}" "" "dot_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}NRM2KERNEL}" "" "nrm2_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}ROTKERNEL}" "" "rot_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}SCALKERNEL}" "" "scal_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}SWAPKERNEL}" "" "swap_k" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}AXPBYKERNEL}" "" "axpby_k" false "" "" false ${float_type})
endforeach ()

# Makefile.L2
GenerateNamedObjects("${KERNELDIR}/gemv_n.S" "DOUBLE")
GenerateNamedObjects("${KERNELDIR}/gemv_t.S" "TRANS")
GenerateCombinationObjects("generic/symv_k.c" "LOWER" "U" "" 1 "" "" 3)
GenerateNamedObjects("generic/ger.c" "" "ger_k" false "" "" "" 3)

# Makefile.L3
foreach (float_type ${FLOAT_TYPES})
  string(SUBSTRING ${float_type} 0 1 float_char)
  GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMKERNEL}" "" "gemm_kernel" false "" "" false ${float_type})
endforeach ()

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
