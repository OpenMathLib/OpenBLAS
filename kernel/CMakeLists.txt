
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

  if (${float_char}GEMMINCOPY)
    GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMINCOPY}" "" "${${float_char}GEMMINCOPYOBJ}" false "" "" false ${float_type})
  endif ()

  if (${float_char}GEMMITCOPY)
    GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMITCOPY}" "" "${${float_char}GEMMITCOPYOBJ}" false "" "" false ${float_type})
  endif ()

  if (${float_char}GEMMONCOPY)
    GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMONCOPY}" "" "${${float_char}GEMMONCOPYOBJ}" false "" "" false ${float_type})
  endif ()

  if (${float_char}GEMMOTCOPY)
    GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMOTCOPY}" "" "${${float_char}GEMMOTCOPYOBJ}" false "" "" false ${float_type})
  endif ()

  GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMM_BETA}" "" "" false "" "" false ${float_type})
  GenerateCombinationObjects("${KERNELDIR}/${${float_char}GEMMKERNEL}" "LEFT;TRANSA" "R;N" "TRMMKERNEL" 2 "trmm_kernel" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMKERNEL}" "UPPER;LN;TRSMKERNEL" "trsm_kernel_LN" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMKERNEL}" "LT;TRSMKERNEL" "trsm_kernel_LT" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMKERNEL}" "UPPER;RN;TRSMKERNEL" "trsm_kernel_RN" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}GEMMKERNEL}" "RT;TRSMKERNEL" "trsm_kernel_RT" false "" "" false ${float_type})

  # These don't use a scheme that is easy to iterate over - the filenames have part of the DEFINE codes in them, for UPPER/TRANS but not for UNIT/OUTER. Also TRANS is not passed in as a define.
  # Could simplify it a bit by pairing up by -UUNIT/-DUNIT.
  GenerateNamedObjects("generic/trmm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "UNIT" "trmm_iunucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "" "trmm_iunncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;UNIT" "trmm_ounucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER" "trmm_ounncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/trmm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER;UNIT" "trmm_ilnucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER" "trmm_ilnncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trmm_olnucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER" "trmm_olnncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/trmm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "UNIT" "trmm_iutucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "" "trmm_iutncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;UNIT" "trmm_outucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER" "trmm_outncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/trmm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER;UNIT" "trmm_iltucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER" "trmm_iltncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trmm_oltucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trmm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER" "trmm_oltncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/trsm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "UNIT" "trsm_iunucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "" "trsm_iunncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;UNIT" "trsm_ounucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_uncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER" "trsm_ounncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/trsm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER;UNIT" "trsm_ilnucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER" "trsm_ilnncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trsm_olnucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_lncopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER" "trsm_olnncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/trsm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "UNIT" "trsm_iutucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "" "trsm_iutncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;UNIT" "trsm_outucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_utcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER" "trsm_outncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/trsm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER;UNIT" "trsm_iltucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "LOWER" "trsm_iltncopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER;UNIT" "trsm_oltucopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/trsm_ltcopy_${${float_char}GEMM_UNROLL_M}.c" "OUTER;LOWER" "trsm_oltncopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/symm_ucopy_${${float_char}GEMM_UNROLL_N}.c" "OUTER" "symm_outcopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/symm_ucopy_${${float_char}GEMM_UNROLL_N}.c" "" "symm_iutcopy" false "" "" false ${float_type})

  GenerateNamedObjects("generic/symm_lcopy_${${float_char}GEMM_UNROLL_N}.c" "LOWER;OUTER" "symm_oltcopy" false "" "" false ${float_type})
  GenerateNamedObjects("generic/symm_lcopy_${${float_char}GEMM_UNROLL_N}.c" "LOWER" "symm_iltcopy" false "" "" false ${float_type})

  if (NOT DEFINED ${float_char}OMATCOPY_CN)
    if (${float_char} STREQUAL "Z" OR ${float_char} STREQUAL "C")
      set(${float_char}OMATCOPY_CN ../arm/zomatcopy_cn.c)
    else ()
      set(${float_char}OMATCOPY_CN ../arm/omatcopy_cn.c)
    endif ()
  endif ()
  if (NOT DEFINED ${float_char}OMATCOPY_RN)
    if (${float_char} STREQUAL "Z" OR ${float_char} STREQUAL "C")
      set(${float_char}OMATCOPY_RN ../arm/zomatcopy_rn.c)
    else ()
      set(${float_char}OMATCOPY_RN ../arm/omatcopy_rn.c)
    endif  ()
  endif ()
  if (NOT DEFINED ${float_char}OMATCOPY_CT)
    if (${float_char} STREQUAL "Z" OR ${float_char} STREQUAL "C")
      set(${float_char}OMATCOPY_CT ../arm/zomatcopy_ct.c)
    else ()
      set(${float_char}OMATCOPY_CT ../arm/omatcopy_ct.c)
    endif ()
  endif ()
  if (NOT DEFINED ${float_char}OMATCOPY_RT)
    if (${float_char} STREQUAL "Z" OR ${float_char} STREQUAL "C")
      set(${float_char}OMATCOPY_RT ../arm/zomatcopy_rt.c)
    else ()
      set(${float_char}OMATCOPY_RT ../arm/omatcopy_rt.c)
    endif ()
  endif ()

  GenerateNamedObjects("${KERNELDIR}/${${float_char}OMATCOPY_CN}" "" "domatcopy_k_cn" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}OMATCOPY_RN}" "ROWM" "domatcopy_k_rn" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}OMATCOPY_CT}" "" "domatcopy_k_ct" false "" "" false ${float_type})
  GenerateNamedObjects("${KERNELDIR}/${${float_char}OMATCOPY_RT}" "ROWM" "domatcopy_k_rt" false "" "" false ${float_type})

endforeach ()

# Makefile.LA
#DBLASOBJS	+= dneg_tcopy$(TSUFFIX).$(SUFFIX) dlaswp_ncopy$(TSUFFIX).$(SUFFIX)

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS
