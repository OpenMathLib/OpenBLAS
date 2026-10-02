/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright 2023, 2025-2026 The OpenBLAS Project.                   */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <stdio.h>
#include <string.h>
#include "common.h"
extern char* gotoblas_corename(void);

#ifdef BUILD_KERNEL
#include "kernelTS.h"
#endif

#undef DEBUG

static void init_parameter(void);

gotoblas_t TABLE_NAME = {
  .core = OPENBLAS_CORETS,
  .dtb_entries = DTB_DEFAULT_ENTRIES,
  .switch_ratio = SWITCH_RATIO,
  .divide_rate = GEMM_DIVIDE_RATE,
  .divide_limit = GEMM_DIVIDE_LIMIT,
  .preferred_size = GEMM_PREFERRED_SIZE,
  .offsetA = GEMM_DEFAULT_OFFSET_A,
  .offsetB = GEMM_DEFAULT_OFFSET_B,
  .align = GEMM_DEFAULT_ALIGN,
#ifdef BUILD_HFLOAT16
  .shgemm_p = 0,
  .shgemm_q = 0,
  .shgemm_r = 0,
  .shgemm_unroll_m = SHGEMM_DEFAULT_UNROLL_M,
  .shgemm_unroll_n = SHGEMM_DEFAULT_UNROLL_N,
#ifdef SHGEMM_DEFAULT_UNROLL_MN
 .shgemm_unroll_mn = SHGEMM_DEFAULT_UNROLL_MN,
#else
 .shgemm_unroll_mn = MAX(SHGEMM_DEFAULT_UNROLL_M, SHGEMM_DEFAULT_UNROLL_N),
#endif
#endif
#ifdef BUILD_BFLOAT16
  .bgemm_p = 0,
  .bgemm_q = 0,
  .bgemm_r = 0,
  .bgemm_unroll_m = BGEMM_DEFAULT_UNROLL_M,
  .bgemm_unroll_n = BGEMM_DEFAULT_UNROLL_N,
#ifdef BGEMM_DEFAULT_UNROLL_MN
 .bgemm_unroll_mn = BGEMM_DEFAULT_UNROLL_MN,
#else
 .bgemm_unroll_mn = MAX(BGEMM_DEFAULT_UNROLL_M, BGEMM_DEFAULT_UNROLL_N),
#endif
  .bgemm_align_k = BGEMM_ALIGN_K,
  .sbgemm_p = 0,
  .sbgemm_q = 0,
  .sbgemm_r = 0,
  .sbgemm_unroll_m = SBGEMM_DEFAULT_UNROLL_M,
  .sbgemm_unroll_n = SBGEMM_DEFAULT_UNROLL_N,
#ifdef SBGEMM_DEFAULT_UNROLL_MN
 .sbgemm_unroll_mn = SBGEMM_DEFAULT_UNROLL_MN,
#else
 .sbgemm_unroll_mn = MAX(SBGEMM_DEFAULT_UNROLL_M, SBGEMM_DEFAULT_UNROLL_N),
#endif
  .sbgemm_align_k = SBGEMM_ALIGN_K,
  .need_amxtile_permission = 0, // need_amxtile_permission
#endif
#if ( BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1) || (BUILD_COMPLEX16==1)
  .sgemm_p = 0,
  .sgemm_q = 0,
  .sgemm_r = 0,
  .sgemm_unroll_m = SGEMM_DEFAULT_UNROLL_M,
  .sgemm_unroll_n = SGEMM_DEFAULT_UNROLL_N,
#ifdef SGEMM_DEFAULT_UNROLL_MN
 .sgemm_unroll_mn = SGEMM_DEFAULT_UNROLL_MN,
#else
 .sgemm_unroll_mn = MAX(SGEMM_DEFAULT_UNROLL_M, SGEMM_DEFAULT_UNROLL_N),
#endif
#endif
#ifdef HAVE_EXCLUSIVE_CACHE
  .exclusive_cache = 1,
#else
  .exclusive_cache = 0,
#endif
#if  (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
  .dgemm_p = 0,
  .dgemm_q = 0,
  .dgemm_r = 0,
  .dgemm_unroll_m = DGEMM_DEFAULT_UNROLL_M,
  .dgemm_unroll_n = DGEMM_DEFAULT_UNROLL_N,
#ifdef DGEMM_DEFAULT_UNROLL_MN
 .dgemm_unroll_mn = DGEMM_DEFAULT_UNROLL_MN,
#else
 .dgemm_unroll_mn = MAX(DGEMM_DEFAULT_UNROLL_M, DGEMM_DEFAULT_UNROLL_N),
#endif
#endif
#ifdef EXPRECISION
  .qgemm_p = 0,
  .qgemm_q = 0,
  .qgemm_r = 0,
  .qgemm_unroll_m = QGEMM_DEFAULT_UNROLL_M,
  .qgemm_unroll_n = QGEMM_DEFAULT_UNROLL_N,
  .qgemm_unroll_mn = MAX(QGEMM_DEFAULT_UNROLL_M, QGEMM_DEFAULT_UNROLL_N),
#endif
#if (BUILD_COMPLEX)
  .cgemm_p = 0,
  .cgemm_q = 0,
  .cgemm_r = 0,
  .cgemm_unroll_m = CGEMM_DEFAULT_UNROLL_M,
  .cgemm_unroll_n = CGEMM_DEFAULT_UNROLL_N,
#ifdef CGEMM_DEFAULT_UNROLL_MN
 .cgemm_unroll_mn = CGEMM_DEFAULT_UNROLL_MN,
#else
 .cgemm_unroll_mn = MAX(CGEMM_DEFAULT_UNROLL_M, CGEMM_DEFAULT_UNROLL_N),
#endif
#endif
#if (BUILD_COMPLEX)
  .cgemm3m_p = 0,
  .cgemm3m_q = 0,
  .cgemm3m_r = 0,
#if (USE_GEMM3M)
#ifdef CGEMM3M_DEFAULT_UNROLL_M
  .cgemm3m_unroll_m = CGEMM3M_DEFAULT_UNROLL_M,
  .cgemm3m_unroll_n = CGEMM3M_DEFAULT_UNROLL_N,
  .cgemm3m_unroll_mn = MAX(CGEMM3M_DEFAULT_UNROLL_M, CGEMM3M_DEFAULT_UNROLL_N),
#else
  .cgemm3m_unroll_m = SGEMM_DEFAULT_UNROLL_M,
  .cgemm3m_unroll_n = SGEMM_DEFAULT_UNROLL_N,
  .cgemm3m_unroll_mn = MAX(SGEMM_DEFAULT_UNROLL_M, SGEMM_DEFAULT_UNROLL_N),
#endif
#else
  .cgemm3m_unroll_m = 0,
  .cgemm3m_unroll_n = 0,
  .cgemm3m_unroll_mn = 0,
#endif
#endif
#if BUILD_COMPLEX16 == 1
  .zgemm_p = 0,
  .zgemm_q = 0,
  .zgemm_r = 0,
  .zgemm_unroll_m = ZGEMM_DEFAULT_UNROLL_M,
  .zgemm_unroll_n = ZGEMM_DEFAULT_UNROLL_N,
#ifdef ZGEMM_DEFAULT_UNROLL_MN
 .zgemm_unroll_mn = ZGEMM_DEFAULT_UNROLL_MN,
#else
 .zgemm_unroll_mn = MAX(ZGEMM_DEFAULT_UNROLL_M, ZGEMM_DEFAULT_UNROLL_N),
#endif
  .zgemm3m_p = 0,
  .zgemm3m_q = 0,
  .zgemm3m_r = 0,
#if (USE_GEMM3M)
#ifdef ZGEMM3M_DEFAULT_UNROLL_M
  .zgemm3m_unroll_m = ZGEMM3M_DEFAULT_UNROLL_M,
  .zgemm3m_unroll_n = ZGEMM3M_DEFAULT_UNROLL_N,
  .zgemm3m_unroll_mn = MAX(ZGEMM3M_DEFAULT_UNROLL_M, ZGEMM3M_DEFAULT_UNROLL_N),
#else
  .zgemm3m_unroll_m = DGEMM_DEFAULT_UNROLL_M,
  .zgemm3m_unroll_n = DGEMM_DEFAULT_UNROLL_N,
  .zgemm3m_unroll_mn = MAX(DGEMM_DEFAULT_UNROLL_M, DGEMM_DEFAULT_UNROLL_N),
#endif
#else
  .zgemm3m_unroll_m = 0,
  .zgemm3m_unroll_n = 0,
  .zgemm3m_unroll_mn = 0,
#endif
#endif
#ifdef EXPRECISION
  .xgemm_p = 0,
  .xgemm_q = 0,
  .xgemm_r = 0,
  .xgemm_unroll_m = XGEMM_DEFAULT_UNROLL_M,
  .xgemm_unroll_n = XGEMM_DEFAULT_UNROLL_N,
  .xgemm_unroll_mn = MAX(XGEMM_DEFAULT_UNROLL_M, XGEMM_DEFAULT_UNROLL_N),
  .xgemm3m_p = 0,
  .xgemm3m_q = 0,
  .xgemm3m_r = 0,
#if (USE_GEMM3M)
  .xgemm3m_unroll_m = QGEMM_DEFAULT_UNROLL_M,
  .xgemm3m_unroll_n = QGEMM_DEFAULT_UNROLL_N,
  .xgemm3m_unroll_mn = MAX(QGEMM_DEFAULT_UNROLL_M, QGEMM_DEFAULT_UNROLL_N),
#else
  .xgemm3m_unroll_m = 0,
  .xgemm3m_unroll_n = 0,
  .xgemm3m_unroll_mn = 0,
#endif
#endif
  .init = init_parameter,
  .snum_opt = SNUMOPT,
  .dnum_opt = DNUMOPT,
  .qnum_opt = QNUMOPT,
};

#if BUILD_HFLOAT16 == 1
const openblas_shgemm_dispatch_t openblas_shgemm_dispatchTS = {
  .shgemm_kernel = shgemm_kernelTS,
  .shgemm_beta = shgemm_betaTS,
#if SHGEMM_DEFAULT_UNROLL_M != SHGEMM_DEFAULT_UNROLL_N
  .shgemm_incopy = shgemm_incopyTS,
  .shgemm_itcopy = shgemm_itcopyTS,
#else
  .shgemm_incopy = shgemm_oncopyTS,
  .shgemm_itcopy = shgemm_otcopyTS,
#endif
  .shgemm_oncopy = shgemm_oncopyTS,
  .shgemm_otcopy = shgemm_otcopyTS,
};
#endif

#if BUILD_HFLOAT16 == 1
const openblas_shgemv_dispatch_t openblas_shgemv_dispatchTS = {
  .shgemv_n = shgemv_nTS,
  .shgemv_t = shgemv_tTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbstobf16_dispatch_t openblas_sbstobf16_dispatchTS = {
  .sbstobf16_k = sbstobf16_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbdtobf16_dispatch_t openblas_sbdtobf16_dispatchTS = {
  .sbdtobf16_k = sbdtobf16_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbf16tos_dispatch_t openblas_sbf16tos_dispatchTS = {
  .sbf16tos_k = sbf16tos_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_dbf16tod_dispatch_t openblas_dbf16tod_dispatchTS = {
  .dbf16tod_k = dbf16tod_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbamax_dispatch_t openblas_sbamax_dispatchTS = {
  .sbamax_k = samax_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbamin_dispatch_t openblas_sbamin_dispatchTS = {
  .sbamin_k = samin_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbmax_dispatch_t openblas_sbmax_dispatchTS = {
  .sbmax_k = smax_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbmin_dispatch_t openblas_sbmin_dispatchTS = {
  .sbmin_k = smin_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_isbamax_dispatch_t openblas_isbamax_dispatchTS = {
  .isbamax_k = isamax_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_isbamin_dispatch_t openblas_isbamin_dispatchTS = {
  .isbamin_k = isamin_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_isbmax_dispatch_t openblas_isbmax_dispatchTS = {
  .isbmax_k = ismax_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_isbmin_dispatch_t openblas_isbmin_dispatchTS = {
  .isbmin_k = ismin_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbnrm2_dispatch_t openblas_sbnrm2_dispatchTS = {
  .sbnrm2_k = snrm2_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbasum_dispatch_t openblas_sbasum_dispatchTS = {
  .sbasum_k = sasum_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbsum_dispatch_t openblas_sbsum_dispatchTS = {
  .sbsum_k = ssum_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbcopy_dispatch_t openblas_sbcopy_dispatchTS = {
  .sbcopy_k = scopy_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbdot_dispatch_t openblas_sbdot_dispatchTS = {
  .sbdot_k = sbdot_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_dsbdot_dispatch_t openblas_dsbdot_dispatchTS = {
  .dsbdot_k = dsdot_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbrot_dispatch_t openblas_sbrot_dispatchTS = {
  .sbrot_k = srot_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbrotm_dispatch_t openblas_sbrotm_dispatchTS = {
  .sbrotm_k = srotm_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_bscal_dispatch_t openblas_bscal_dispatchTS = {
  .bscal_k = bscal_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbaxpy_dispatch_t openblas_sbaxpy_dispatchTS = {
  .sbaxpy_k = saxpy_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbscal_dispatch_t openblas_sbscal_dispatchTS = {
  .sbscal_k = sscal_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbswap_dispatch_t openblas_sbswap_dispatchTS = {
  .sbswap_k = sswap_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_bgemv_dispatch_t openblas_bgemv_dispatchTS = {
  .bgemv_n = bgemv_nTS,
  .bgemv_t = bgemv_tTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbgemv_dispatch_t openblas_sbgemv_dispatchTS = {
  .sbgemv_n = sbgemv_nTS,
  .sbgemv_t = sbgemv_tTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbger_dispatch_t openblas_sbger_dispatchTS = {
  .sbger_k = sger_kTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbsymv_dispatch_t openblas_sbsymv_dispatchTS = {
  .sbsymv_L = ssymv_LTS,
  .sbsymv_U = ssymv_UTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_bgemm_dispatch_t openblas_bgemm_dispatchTS = {
  .bgemm_kernel = bgemm_kernelTS,
  .bgemm_beta = bgemm_betaTS,
#if BGEMM_DEFAULT_UNROLL_M != BGEMM_DEFAULT_UNROLL_N
  .bgemm_incopy = bgemm_incopyTS,
  .bgemm_itcopy = bgemm_itcopyTS,
#else
  .bgemm_incopy = bgemm_oncopyTS,
  .bgemm_itcopy = bgemm_otcopyTS,
#endif
  .bgemm_oncopy = bgemm_oncopyTS,
  .bgemm_otcopy = bgemm_otcopyTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbgemm_dispatch_t openblas_sbgemm_dispatchTS = {
  .sbgemm_kernel = sbgemm_kernelTS,
  .sbgemm_beta = sbgemm_betaTS,
#if SBGEMM_DEFAULT_UNROLL_M != SBGEMM_DEFAULT_UNROLL_N
  .sbgemm_incopy = sbgemm_incopyTS,
  .sbgemm_itcopy = sbgemm_itcopyTS,
#else
  .sbgemm_incopy = sbgemm_oncopyTS,
  .sbgemm_itcopy = sbgemm_otcopyTS,
#endif
  .sbgemm_oncopy = sbgemm_oncopyTS,
  .sbgemm_otcopy = sbgemm_otcopyTS,
#ifdef SMALL_MATRIX_OPT
  .sbgemm_small_matrix_permit = sbgemm_small_matrix_permitTS,
  .sbgemm_small_kernel_nn = sbgemm_small_kernel_nnTS,
  .sbgemm_small_kernel_nt = sbgemm_small_kernel_ntTS,
  .sbgemm_small_kernel_tn = sbgemm_small_kernel_tnTS,
  .sbgemm_small_kernel_tt = sbgemm_small_kernel_ttTS,
  .sbgemm_small_kernel_b0_nn = sbgemm_small_kernel_b0_nnTS,
  .sbgemm_small_kernel_b0_nt = sbgemm_small_kernel_b0_ntTS,
  .sbgemm_small_kernel_b0_tn = sbgemm_small_kernel_b0_tnTS,
  .sbgemm_small_kernel_b0_tt = sbgemm_small_kernel_b0_ttTS,
#endif
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbtrsm_dispatch_t openblas_sbtrsm_dispatchTS = {
  .sbtrsm_kernel_LN = strsm_kernel_LNTS,
  .sbtrsm_kernel_LT = strsm_kernel_LTTS,
  .sbtrsm_kernel_RN = strsm_kernel_RNTS,
  .sbtrsm_kernel_RT = strsm_kernel_RTTS,
#if SGEMM_DEFAULT_UNROLL_M != SGEMM_DEFAULT_UNROLL_N
  .sbtrsm_iunucopy = strsm_iunucopyTS,
  .sbtrsm_iunncopy = strsm_iunncopyTS,
  .sbtrsm_iutucopy = strsm_iutucopyTS,
  .sbtrsm_iutncopy = strsm_iutncopyTS,
  .sbtrsm_ilnucopy = strsm_ilnucopyTS,
  .sbtrsm_ilnncopy = strsm_ilnncopyTS,
  .sbtrsm_iltucopy = strsm_iltucopyTS,
  .sbtrsm_iltncopy = strsm_iltncopyTS,
#else
  .sbtrsm_iunucopy = strsm_ounucopyTS,
  .sbtrsm_iunncopy = strsm_ounncopyTS,
  .sbtrsm_iutucopy = strsm_outucopyTS,
  .sbtrsm_iutncopy = strsm_outncopyTS,
  .sbtrsm_ilnucopy = strsm_olnucopyTS,
  .sbtrsm_ilnncopy = strsm_olnncopyTS,
  .sbtrsm_iltucopy = strsm_oltucopyTS,
  .sbtrsm_iltncopy = strsm_oltncopyTS,
#endif
  .sbtrsm_ounucopy = strsm_ounucopyTS,
  .sbtrsm_ounncopy = strsm_ounncopyTS,
  .sbtrsm_outucopy = strsm_outucopyTS,
  .sbtrsm_outncopy = strsm_outncopyTS,
  .sbtrsm_olnucopy = strsm_olnucopyTS,
  .sbtrsm_olnncopy = strsm_olnncopyTS,
  .sbtrsm_oltucopy = strsm_oltucopyTS,
  .sbtrsm_oltncopy = strsm_oltncopyTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbtrmm_dispatch_t openblas_sbtrmm_dispatchTS = {
  .sbtrmm_kernel_RN = strmm_kernel_RNTS,
  .sbtrmm_kernel_RT = strmm_kernel_RTTS,
  .sbtrmm_kernel_LN = strmm_kernel_LNTS,
  .sbtrmm_kernel_LT = strmm_kernel_LTTS,
#if SGEMM_DEFAULT_UNROLL_M != SGEMM_DEFAULT_UNROLL_N
  .sbtrmm_iunucopy = strmm_iunucopyTS,
  .sbtrmm_iunncopy = strmm_iunncopyTS,
  .sbtrmm_iutucopy = strmm_iutucopyTS,
  .sbtrmm_iutncopy = strmm_iutncopyTS,
  .sbtrmm_ilnucopy = strmm_ilnucopyTS,
  .sbtrmm_ilnncopy = strmm_ilnncopyTS,
  .sbtrmm_iltucopy = strmm_iltucopyTS,
  .sbtrmm_iltncopy = strmm_iltncopyTS,
#else
  .sbtrmm_iunucopy = strmm_ounucopyTS,
  .sbtrmm_iunncopy = strmm_ounncopyTS,
  .sbtrmm_iutucopy = strmm_outucopyTS,
  .sbtrmm_iutncopy = strmm_outncopyTS,
  .sbtrmm_ilnucopy = strmm_olnucopyTS,
  .sbtrmm_ilnncopy = strmm_olnncopyTS,
  .sbtrmm_iltucopy = strmm_oltucopyTS,
  .sbtrmm_iltncopy = strmm_oltncopyTS,
#endif
  .sbtrmm_ounucopy = strmm_ounucopyTS,
  .sbtrmm_ounncopy = strmm_ounncopyTS,
  .sbtrmm_outucopy = strmm_outucopyTS,
  .sbtrmm_outncopy = strmm_outncopyTS,
  .sbtrmm_olnucopy = strmm_olnucopyTS,
  .sbtrmm_olnncopy = strmm_olnncopyTS,
  .sbtrmm_oltucopy = strmm_oltucopyTS,
  .sbtrmm_oltncopy = strmm_oltncopyTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbsymm_dispatch_t openblas_sbsymm_dispatchTS = {
#if SGEMM_DEFAULT_UNROLL_M != SGEMM_DEFAULT_UNROLL_N
  .sbsymm_iutcopy = ssymm_iutcopyTS,
  .sbsymm_iltcopy = ssymm_iltcopyTS,
#else
  .sbsymm_iutcopy = ssymm_outcopyTS,
  .sbsymm_iltcopy = ssymm_oltcopyTS,
#endif
  .sbsymm_outcopy = ssymm_outcopyTS,
  .sbsymm_oltcopy = ssymm_oltcopyTS,
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sbneg_dispatch_t openblas_sbneg_dispatchTS = {
#ifndef NO_LAPACK
  .sbneg_tcopy = sneg_tcopyTS,
#else
  .sbneg_tcopy = NULL,
#endif
};
#endif

#if BUILD_BFLOAT16 == 1
const openblas_sblaswp_dispatch_t openblas_sblaswp_dispatchTS = {
#ifndef NO_LAPACK
  .sblaswp_ncopy = slaswp_ncopyTS,
#else
  .sblaswp_ncopy = NULL,
#endif
};
#endif

#if (BUILD_SINGLE == 1) || (BUILD_COMPLEX == 1)
const openblas_samax_dispatch_t openblas_samax_dispatchTS = {
  .samax_k = samax_kTS,
};
#endif

#if (BUILD_SINGLE == 1) || (BUILD_COMPLEX == 1)
const openblas_samin_dispatch_t openblas_samin_dispatchTS = {
  .samin_k = samin_kTS,
};
#endif

#if (BUILD_SINGLE == 1) || (BUILD_COMPLEX == 1)
const openblas_smax_dispatch_t openblas_smax_dispatchTS = {
  .smax_k = smax_kTS,
};
#endif

#if (BUILD_SINGLE == 1) || (BUILD_COMPLEX == 1)
const openblas_smin_dispatch_t openblas_smin_dispatchTS = {
  .smin_k = smin_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE ==1) || (BUILD_COMPLEX==1)
const openblas_isamax_dispatch_t openblas_isamax_dispatchTS = {
  .isamax_k = isamax_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
const openblas_isamin_dispatch_t openblas_isamin_dispatchTS = {
  .isamin_k = isamin_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
const openblas_ismax_dispatch_t openblas_ismax_dispatchTS = {
  .ismax_k = ismax_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
const openblas_ismin_dispatch_t openblas_ismin_dispatchTS = {
  .ismin_k = ismin_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
const openblas_snrm2_dispatch_t openblas_snrm2_dispatchTS = {
  .snrm2_k = snrm2_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
const openblas_sasum_dispatch_t openblas_sasum_dispatchTS = {
  .sasum_k = sasum_kTS,
};
#endif

#if (BUILD_SINGLE==1)
const openblas_ssum_dispatch_t openblas_ssum_dispatchTS = {
  .ssum_k = ssum_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_scopy_dispatch_t openblas_scopy_dispatchTS = {
  .scopy_k = scopy_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_sdot_dispatch_t openblas_sdot_dispatchTS = {
  .sdot_k = sdot_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_srot_dispatch_t openblas_srot_dispatchTS = {
  .srot_k = srot_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_srotm_dispatch_t openblas_srotm_dispatchTS = {
  .srotm_k = srotm_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_saxpy_dispatch_t openblas_saxpy_dispatchTS = {
  .saxpy_k = saxpy_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1) || (BUILD_COMPLEX16==1)
const openblas_sscal_dispatch_t openblas_sscal_dispatchTS = {
  .sscal_k = sscal_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_sswap_dispatch_t openblas_sswap_dispatchTS = {
  .sswap_k = sswap_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_sgemv_dispatch_t openblas_sgemv_dispatchTS = {
  .sgemv_n = sgemv_nTS,
  .sgemv_t = sgemv_tTS,
};
#endif

#if (BUILD_SINGLE==1)
const openblas_sger_dispatch_t openblas_sger_dispatchTS = {
  .sger_k = sger_kTS,
};
#endif

#if (BUILD_SINGLE==1)
const openblas_ssymv_dispatch_t openblas_ssymv_dispatchTS = {
  .ssymv_L = ssymv_LTS,
  .ssymv_U = ssymv_UTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_sgemm_dispatch_t openblas_sgemm_dispatchTS = {
#ifdef ARCH_X86_64
  .sgemm_direct = sgemm_directTS,
  .sgemm_direct_performant = sgemm_direct_performantTS,
#endif
#ifdef ARCH_ARM64
  .sgemm_direct = sgemm_directTS,
  .sgemm_direct_performant = sgemm_direct_performantTS,
  .sgemm_direct_alpha_beta = sgemm_direct_alpha_betaTS,
#ifdef HAVE_SME
  .sme_sgemm_kernel = sme_sgemm_kernelTS,
#else
  .sme_sgemm_kernel = NULL,
#endif
#endif
  .sgemm_kernel = sgemm_kernelTS,
  .sgemm_beta = sgemm_betaTS,
#if SGEMM_DEFAULT_UNROLL_M != SGEMM_DEFAULT_UNROLL_N
  .sgemm_incopy = sgemm_incopyTS,
  .sgemm_itcopy = sgemm_itcopyTS,
#else
  .sgemm_incopy = sgemm_oncopyTS,
  .sgemm_itcopy = sgemm_otcopyTS,
#endif
  .sgemm_oncopy = sgemm_oncopyTS,
  .sgemm_otcopy = sgemm_otcopyTS,
#ifdef SMALL_MATRIX_OPT
  .sgemm_small_matrix_permit = sgemm_small_matrix_permitTS,
  .sgemm_small_kernel_nn = sgemm_small_kernel_nnTS,
  .sgemm_small_kernel_nt = sgemm_small_kernel_ntTS,
  .sgemm_small_kernel_tn = sgemm_small_kernel_tnTS,
  .sgemm_small_kernel_tt = sgemm_small_kernel_ttTS,
  .sgemm_small_kernel_b0_nn = sgemm_small_kernel_b0_nnTS,
  .sgemm_small_kernel_b0_nt = sgemm_small_kernel_b0_ntTS,
  .sgemm_small_kernel_b0_tn = sgemm_small_kernel_b0_tnTS,
  .sgemm_small_kernel_b0_tt = sgemm_small_kernel_b0_ttTS,
#endif
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_ssymm_dispatch_t openblas_ssymm_dispatchTS = {
#ifdef ARCH_ARM64
  .ssymm_direct_alpha_betaLU = ssymm_direct_alpha_betaLUTS,
  .ssymm_direct_alpha_betaLL = ssymm_direct_alpha_betaLLTS,
#endif
  .ssymm_kernel = ssymm_kernelTS,
#if (BUILD_SINGLE==1)
  .ssymm_incopy = ssymm_incopyTS,
  .ssymm_itcopy = ssymm_itcopyTS,
#if SGEMM_DEFAULT_UNROLL_M != SGEMM_DEFAULT_UNROLL_N
  .ssymm_iutcopy = ssymm_iutcopyTS,
  .ssymm_iltcopy = ssymm_iltcopyTS,
#else
  .ssymm_iutcopy = ssymm_outcopyTS,
  .ssymm_iltcopy = ssymm_oltcopyTS,
#endif
  .ssymm_outcopy = ssymm_outcopyTS,
  .ssymm_oltcopy = ssymm_oltcopyTS,
#endif
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_strmm_dispatch_t openblas_strmm_dispatchTS = {
#ifdef ARCH_ARM64
  .strmm_direct_LNUN = strmm_direct_LNUNTS,
  .strmm_direct_LNLN = strmm_direct_LNLNTS,
  .strmm_direct_LTUN = strmm_direct_LTUNTS,
  .strmm_direct_LTLN = strmm_direct_LTLNTS,
#endif
  .strmm_gemm_kernel = strmm_gemm_kernelTS,
#if (BUILD_SINGLE==1)
  .strmm_kernel_RN = strmm_kernel_RNTS,
  .strmm_kernel_RT = strmm_kernel_RTTS,
  .strmm_kernel_LN = strmm_kernel_LNTS,
  .strmm_kernel_LT = strmm_kernel_LTTS,
#if SGEMM_DEFAULT_UNROLL_M != SGEMM_DEFAULT_UNROLL_N
  .strmm_iunucopy = strmm_iunucopyTS,
  .strmm_iunncopy = strmm_iunncopyTS,
  .strmm_iutucopy = strmm_iutucopyTS,
  .strmm_iutncopy = strmm_iutncopyTS,
  .strmm_ilnucopy = strmm_ilnucopyTS,
  .strmm_ilnncopy = strmm_ilnncopyTS,
  .strmm_iltucopy = strmm_iltucopyTS,
  .strmm_iltncopy = strmm_iltncopyTS,
#else
  .strmm_iunucopy = strmm_ounucopyTS,
  .strmm_iunncopy = strmm_ounncopyTS,
  .strmm_iutucopy = strmm_outucopyTS,
  .strmm_iutncopy = strmm_outncopyTS,
  .strmm_ilnucopy = strmm_olnucopyTS,
  .strmm_ilnncopy = strmm_olnncopyTS,
  .strmm_iltucopy = strmm_oltucopyTS,
  .strmm_iltncopy = strmm_oltncopyTS,
#endif
  .strmm_incopy = strmm_incopyTS,
  .strmm_itcopy = strmm_itcopyTS,
  .strmm_ounucopy = strmm_ounucopyTS,
  .strmm_ounncopy = strmm_ounncopyTS,
  .strmm_outucopy = strmm_outucopyTS,
  .strmm_outncopy = strmm_outncopyTS,
  .strmm_olnucopy = strmm_olnucopyTS,
  .strmm_olnncopy = strmm_olnncopyTS,
  .strmm_oltucopy = strmm_oltucopyTS,
  .strmm_oltncopy = strmm_oltncopyTS,
#endif
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
#ifdef ARCH_ARM64
const openblas_ssyrk_dispatch_t openblas_ssyrk_dispatchTS = {
  .ssyrk_direct_alpha_betaUN = ssyrk_direct_alpha_betaUNTS,
  .ssyrk_direct_alpha_betaUT = ssyrk_direct_alpha_betaUTTS,
  .ssyrk_direct_alpha_betaLN = ssyrk_direct_alpha_betaLNTS,
  .ssyrk_direct_alpha_betaLT = ssyrk_direct_alpha_betaLTTS,
};
#endif
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
#ifdef ARCH_ARM64
const openblas_ssyr2k_dispatch_t openblas_ssyr2k_dispatchTS = {
  .ssyr2k_direct_alpha_betaUN = ssyr2k_direct_alpha_betaUNTS,
  .ssyr2k_direct_alpha_betaUT = ssyr2k_direct_alpha_betaUTTS,
  .ssyr2k_direct_alpha_betaLN = ssyr2k_direct_alpha_betaLNTS,
  .ssyr2k_direct_alpha_betaLT = ssyr2k_direct_alpha_betaLTTS,
};
#endif
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
const openblas_strsm_dispatch_t openblas_strsm_dispatchTS = {
  .strsm_kernel_LN = strsm_kernel_LNTS,
  .strsm_kernel_LT = strsm_kernel_LTTS,
  .strsm_kernel_RN = strsm_kernel_RNTS,
  .strsm_kernel_RT = strsm_kernel_RTTS,
#if SGEMM_DEFAULT_UNROLL_M != SGEMM_DEFAULT_UNROLL_N
  .strsm_iunucopy = strsm_iunucopyTS,
  .strsm_iunncopy = strsm_iunncopyTS,
  .strsm_iutucopy = strsm_iutucopyTS,
  .strsm_iutncopy = strsm_iutncopyTS,
  .strsm_ilnucopy = strsm_ilnucopyTS,
  .strsm_ilnncopy = strsm_ilnncopyTS,
  .strsm_iltucopy = strsm_iltucopyTS,
  .strsm_iltncopy = strsm_iltncopyTS,
#else
  .strsm_iunucopy = strsm_ounucopyTS,
  .strsm_iunncopy = strsm_ounncopyTS,
  .strsm_iutucopy = strsm_outucopyTS,
  .strsm_iutncopy = strsm_outncopyTS,
  .strsm_ilnucopy = strsm_olnucopyTS,
  .strsm_ilnncopy = strsm_olnncopyTS,
  .strsm_iltucopy = strsm_oltucopyTS,
  .strsm_iltncopy = strsm_oltncopyTS,
#endif
  .strsm_ounucopy = strsm_ounucopyTS,
  .strsm_ounncopy = strsm_ounncopyTS,
  .strsm_outucopy = strsm_outucopyTS,
  .strsm_outncopy = strsm_outncopyTS,
  .strsm_olnucopy = strsm_olnucopyTS,
  .strsm_olnncopy = strsm_olnncopyTS,
  .strsm_oltucopy = strsm_oltucopyTS,
  .strsm_oltncopy = strsm_oltncopyTS,
};
#endif

#if (BUILD_SINGLE==1)
const openblas_sneg_dispatch_t openblas_sneg_dispatchTS = {
#ifndef NO_LAPACK
  .sneg_tcopy = sneg_tcopyTS,
#else
  .sneg_tcopy = NULL,
#endif
};
#endif

#if (BUILD_SINGLE==1)
const openblas_slaswp_dispatch_t openblas_slaswp_dispatchTS = {
#ifndef NO_LAPACK
  .slaswp_ncopy = slaswp_ncopyTS,
#else
  .slaswp_ncopy = NULL,
#endif
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_damax_dispatch_t openblas_damax_dispatchTS = {
  .damax_k = damax_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_damin_dispatch_t openblas_damin_dispatchTS = {
  .damin_k = damin_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dmax_dispatch_t openblas_dmax_dispatchTS = {
  .dmax_k = dmax_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dmin_dispatch_t openblas_dmin_dispatchTS = {
  .dmin_k = dmin_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_idamax_dispatch_t openblas_idamax_dispatchTS = {
  .idamax_k = idamax_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_idamin_dispatch_t openblas_idamin_dispatchTS = {
  .idamin_k = idamin_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_idmax_dispatch_t openblas_idmax_dispatchTS = {
  .idmax_k = idmax_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_idmin_dispatch_t openblas_idmin_dispatchTS = {
  .idmin_k = idmin_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dnrm2_dispatch_t openblas_dnrm2_dispatchTS = {
  .dnrm2_k = dnrm2_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dasum_dispatch_t openblas_dasum_dispatchTS = {
  .dasum_k = dasum_kTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dsum_dispatch_t openblas_dsum_dispatchTS = {
  .dsum_k = dsum_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dcopy_dispatch_t openblas_dcopy_dispatchTS = {
  .dcopy_k = dcopy_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_ddot_dispatch_t openblas_ddot_dispatchTS = {
  .ddot_k = ddot_kTS,
};
#endif

#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1)
const openblas_dsdot_dispatch_t openblas_dsdot_dispatchTS = {
  .dsdot_k = dsdot_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_drot_dispatch_t openblas_drot_dispatchTS = {
  .drot_k = drot_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_drotm_dispatch_t openblas_drotm_dispatchTS = {
  .drotm_k = drotm_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_daxpy_dispatch_t openblas_daxpy_dispatchTS = {
  .daxpy_k = daxpy_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dscal_dispatch_t openblas_dscal_dispatchTS = {
  .dscal_k = dscal_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dswap_dispatch_t openblas_dswap_dispatchTS = {
  .dswap_k = dswap_kTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dgemv_dispatch_t openblas_dgemv_dispatchTS = {
  .dgemv_n = dgemv_nTS,
  .dgemv_t = dgemv_tTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dger_dispatch_t openblas_dger_dispatchTS = {
  .dger_k = dger_kTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dsymv_dispatch_t openblas_dsymv_dispatchTS = {
  .dsymv_L = dsymv_LTS,
  .dsymv_U = dsymv_UTS,
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dgemm_dispatch_t openblas_dgemm_dispatchTS = {
#ifdef ARCH_ARM64
#ifdef HAVE_SME
  .sme_dgemm_kernel = sme_dgemm_kernelTS,
#else
  .sme_dgemm_kernel = NULL,
#endif
#endif
  .dgemm_kernel = dgemm_kernelTS,
  .dgemm_beta = dgemm_betaTS,
#if DGEMM_DEFAULT_UNROLL_M != DGEMM_DEFAULT_UNROLL_N
  .dgemm_incopy = dgemm_incopyTS,
  .dgemm_itcopy = dgemm_itcopyTS,
#else
  .dgemm_incopy = dgemm_oncopyTS,
  .dgemm_itcopy = dgemm_otcopyTS,
#endif
  .dgemm_oncopy = dgemm_oncopyTS,
  .dgemm_otcopy = dgemm_otcopyTS,
#ifdef SMALL_MATRIX_OPT
  .dgemm_small_matrix_permit = dgemm_small_matrix_permitTS,
  .dgemm_small_kernel_nn = dgemm_small_kernel_nnTS,
  .dgemm_small_kernel_nt = dgemm_small_kernel_ntTS,
  .dgemm_small_kernel_tn = dgemm_small_kernel_tnTS,
  .dgemm_small_kernel_tt = dgemm_small_kernel_ttTS,
  .dgemm_small_kernel_b0_nn = dgemm_small_kernel_b0_nnTS,
  .dgemm_small_kernel_b0_nt = dgemm_small_kernel_b0_ntTS,
  .dgemm_small_kernel_b0_tn = dgemm_small_kernel_b0_tnTS,
  .dgemm_small_kernel_b0_tt = dgemm_small_kernel_b0_ttTS,
#endif
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dsymm_dispatch_t openblas_dsymm_dispatchTS = {
  .dsymm_kernel = dsymm_kernelTS,
#if  (BUILD_DOUBLE==1)
  .dsymm_incopy = dsymm_incopyTS,
  .dsymm_itcopy = dsymm_itcopyTS,
#if DGEMM_DEFAULT_UNROLL_M != DGEMM_DEFAULT_UNROLL_N
  .dsymm_iutcopy = dsymm_iutcopyTS,
  .dsymm_iltcopy = dsymm_iltcopyTS,
#else
  .dsymm_iutcopy = dsymm_outcopyTS,
  .dsymm_iltcopy = dsymm_oltcopyTS,
#endif
  .dsymm_outcopy = dsymm_outcopyTS,
  .dsymm_oltcopy = dsymm_oltcopyTS,
#endif
};
#endif

#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
const openblas_dtrmm_dispatch_t openblas_dtrmm_dispatchTS = {
  .dtrmm_gemm_kernel = dtrmm_gemm_kernelTS,
#if  (BUILD_DOUBLE==1)
  .dtrmm_kernel_RN = dtrmm_kernel_RNTS,
  .dtrmm_kernel_RT = dtrmm_kernel_RTTS,
  .dtrmm_kernel_LN = dtrmm_kernel_LNTS,
  .dtrmm_kernel_LT = dtrmm_kernel_LTTS,
#if DGEMM_DEFAULT_UNROLL_M != DGEMM_DEFAULT_UNROLL_N
  .dtrmm_iunucopy = dtrmm_iunucopyTS,
  .dtrmm_iunncopy = dtrmm_iunncopyTS,
  .dtrmm_iutucopy = dtrmm_iutucopyTS,
  .dtrmm_iutncopy = dtrmm_iutncopyTS,
  .dtrmm_ilnucopy = dtrmm_ilnucopyTS,
  .dtrmm_ilnncopy = dtrmm_ilnncopyTS,
  .dtrmm_iltucopy = dtrmm_iltucopyTS,
  .dtrmm_iltncopy = dtrmm_iltncopyTS,
#else
  .dtrmm_iunucopy = dtrmm_ounucopyTS,
  .dtrmm_iunncopy = dtrmm_ounncopyTS,
  .dtrmm_iutucopy = dtrmm_outucopyTS,
  .dtrmm_iutncopy = dtrmm_outncopyTS,
  .dtrmm_ilnucopy = dtrmm_olnucopyTS,
  .dtrmm_ilnncopy = dtrmm_olnncopyTS,
  .dtrmm_iltucopy = dtrmm_oltucopyTS,
  .dtrmm_iltncopy = dtrmm_oltncopyTS,
#endif
  .dtrmm_incopy = dtrmm_incopyTS,
  .dtrmm_itcopy = dtrmm_itcopyTS,
  .dtrmm_ounucopy = dtrmm_ounucopyTS,
  .dtrmm_ounncopy = dtrmm_ounncopyTS,
  .dtrmm_outucopy = dtrmm_outucopyTS,
  .dtrmm_outncopy = dtrmm_outncopyTS,
  .dtrmm_olnucopy = dtrmm_olnucopyTS,
  .dtrmm_olnncopy = dtrmm_olnncopyTS,
  .dtrmm_oltucopy = dtrmm_oltucopyTS,
  .dtrmm_oltncopy = dtrmm_oltncopyTS,
#endif
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dtrsm_dispatch_t openblas_dtrsm_dispatchTS = {
  .dtrsm_kernel_LN = dtrsm_kernel_LNTS,
  .dtrsm_kernel_LT = dtrsm_kernel_LTTS,
  .dtrsm_kernel_RN = dtrsm_kernel_RNTS,
  .dtrsm_kernel_RT = dtrsm_kernel_RTTS,
#if DGEMM_DEFAULT_UNROLL_M != DGEMM_DEFAULT_UNROLL_N
  .dtrsm_iunucopy = dtrsm_iunucopyTS,
  .dtrsm_iunncopy = dtrsm_iunncopyTS,
  .dtrsm_iutucopy = dtrsm_iutucopyTS,
  .dtrsm_iutncopy = dtrsm_iutncopyTS,
  .dtrsm_ilnucopy = dtrsm_ilnucopyTS,
  .dtrsm_ilnncopy = dtrsm_ilnncopyTS,
  .dtrsm_iltucopy = dtrsm_iltucopyTS,
  .dtrsm_iltncopy = dtrsm_iltncopyTS,
#else
  .dtrsm_iunucopy = dtrsm_ounucopyTS,
  .dtrsm_iunncopy = dtrsm_ounncopyTS,
  .dtrsm_iutucopy = dtrsm_outucopyTS,
  .dtrsm_iutncopy = dtrsm_outncopyTS,
  .dtrsm_ilnucopy = dtrsm_olnucopyTS,
  .dtrsm_ilnncopy = dtrsm_olnncopyTS,
  .dtrsm_iltucopy = dtrsm_oltucopyTS,
  .dtrsm_iltncopy = dtrsm_oltncopyTS,
#endif
  .dtrsm_ounucopy = dtrsm_ounucopyTS,
  .dtrsm_ounncopy = dtrsm_ounncopyTS,
  .dtrsm_outucopy = dtrsm_outucopyTS,
  .dtrsm_outncopy = dtrsm_outncopyTS,
  .dtrsm_olnucopy = dtrsm_olnucopyTS,
  .dtrsm_olnncopy = dtrsm_olnncopyTS,
  .dtrsm_oltucopy = dtrsm_oltucopyTS,
  .dtrsm_oltncopy = dtrsm_oltncopyTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dneg_dispatch_t openblas_dneg_dispatchTS = {
#ifndef NO_LAPACK
  .dneg_tcopy = dneg_tcopyTS,
#else
  .dneg_tcopy = NULL,
#endif
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dlaswp_dispatch_t openblas_dlaswp_dispatchTS = {
#ifndef NO_LAPACK
  .dlaswp_ncopy = dlaswp_ncopyTS,
#else
  .dlaswp_ncopy = NULL,
#endif
};
#endif

#ifdef EXPRECISION
const openblas_qamax_dispatch_t openblas_qamax_dispatchTS = {
  .qamax_k = qamax_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qamin_dispatch_t openblas_qamin_dispatchTS = {
  .qamin_k = qamin_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qmax_dispatch_t openblas_qmax_dispatchTS = {
  .qmax_k = qmax_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qmin_dispatch_t openblas_qmin_dispatchTS = {
  .qmin_k = qmin_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_iqamax_dispatch_t openblas_iqamax_dispatchTS = {
  .iqamax_k = iqamax_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_iqamin_dispatch_t openblas_iqamin_dispatchTS = {
  .iqamin_k = iqamin_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_iqmax_dispatch_t openblas_iqmax_dispatchTS = {
  .iqmax_k = iqmax_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_iqmin_dispatch_t openblas_iqmin_dispatchTS = {
  .iqmin_k = iqmin_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qnrm2_dispatch_t openblas_qnrm2_dispatchTS = {
  .qnrm2_k = qnrm2_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qasum_dispatch_t openblas_qasum_dispatchTS = {
  .qasum_k = qasum_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qsum_dispatch_t openblas_qsum_dispatchTS = {
  .qsum_k = qsum_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qcopy_dispatch_t openblas_qcopy_dispatchTS = {
  .qcopy_k = qcopy_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qdot_dispatch_t openblas_qdot_dispatchTS = {
  .qdot_k = qdot_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qrot_dispatch_t openblas_qrot_dispatchTS = {
  .qrot_k = qrot_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qrotm_dispatch_t openblas_qrotm_dispatchTS = {
  .qrotm_k = qrotm_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qaxpy_dispatch_t openblas_qaxpy_dispatchTS = {
  .qaxpy_k = qaxpy_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qscal_dispatch_t openblas_qscal_dispatchTS = {
  .qscal_k = qscal_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qswap_dispatch_t openblas_qswap_dispatchTS = {
  .qswap_k = qswap_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qgemv_dispatch_t openblas_qgemv_dispatchTS = {
  .qgemv_n = qgemv_nTS,
  .qgemv_t = qgemv_tTS,
};
#endif

#ifdef EXPRECISION
const openblas_qger_dispatch_t openblas_qger_dispatchTS = {
  .qger_k = qger_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_qsymv_dispatch_t openblas_qsymv_dispatchTS = {
  .qsymv_L = qsymv_LTS,
  .qsymv_U = qsymv_UTS,
};
#endif

#ifdef EXPRECISION
const openblas_qgemm_dispatch_t openblas_qgemm_dispatchTS = {
  .qgemm_kernel = qgemm_kernelTS,
  .qgemm_beta = qgemm_betaTS,
#if QGEMM_DEFAULT_UNROLL_M != QGEMM_DEFAULT_UNROLL_N
  .qgemm_incopy = qgemm_incopyTS,
  .qgemm_itcopy = qgemm_itcopyTS,
#else
  .qgemm_incopy = qgemm_oncopyTS,
  .qgemm_itcopy = qgemm_otcopyTS,
#endif
  .qgemm_oncopy = qgemm_oncopyTS,
  .qgemm_otcopy = qgemm_otcopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_qtrsm_dispatch_t openblas_qtrsm_dispatchTS = {
  .qtrsm_kernel_LN = qtrsm_kernel_LNTS,
  .qtrsm_kernel_LT = qtrsm_kernel_LTTS,
  .qtrsm_kernel_RN = qtrsm_kernel_RNTS,
  .qtrsm_kernel_RT = qtrsm_kernel_RTTS,
#if QGEMM_DEFAULT_UNROLL_M != QGEMM_DEFAULT_UNROLL_N
  .qtrsm_iunucopy = qtrsm_iunucopyTS,
  .qtrsm_iunncopy = qtrsm_iunncopyTS,
  .qtrsm_iutucopy = qtrsm_iutucopyTS,
  .qtrsm_iutncopy = qtrsm_iutncopyTS,
  .qtrsm_ilnucopy = qtrsm_ilnucopyTS,
  .qtrsm_ilnncopy = qtrsm_ilnncopyTS,
  .qtrsm_iltucopy = qtrsm_iltucopyTS,
  .qtrsm_iltncopy = qtrsm_iltncopyTS,
#else
  .qtrsm_iunucopy = qtrsm_ounucopyTS,
  .qtrsm_iunncopy = qtrsm_ounncopyTS,
  .qtrsm_iutucopy = qtrsm_outucopyTS,
  .qtrsm_iutncopy = qtrsm_outncopyTS,
  .qtrsm_ilnucopy = qtrsm_olnucopyTS,
  .qtrsm_ilnncopy = qtrsm_olnncopyTS,
  .qtrsm_iltucopy = qtrsm_oltucopyTS,
  .qtrsm_iltncopy = qtrsm_oltncopyTS,
#endif
  .qtrsm_ounucopy = qtrsm_ounucopyTS,
  .qtrsm_ounncopy = qtrsm_ounncopyTS,
  .qtrsm_outucopy = qtrsm_outucopyTS,
  .qtrsm_outncopy = qtrsm_outncopyTS,
  .qtrsm_olnucopy = qtrsm_olnucopyTS,
  .qtrsm_olnncopy = qtrsm_olnncopyTS,
  .qtrsm_oltucopy = qtrsm_oltucopyTS,
  .qtrsm_oltncopy = qtrsm_oltncopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_qtrmm_dispatch_t openblas_qtrmm_dispatchTS = {
  .qtrmm_kernel_RN = qtrmm_kernel_RNTS,
  .qtrmm_kernel_RT = qtrmm_kernel_RTTS,
  .qtrmm_kernel_LN = qtrmm_kernel_LNTS,
  .qtrmm_kernel_LT = qtrmm_kernel_LTTS,
#if QGEMM_DEFAULT_UNROLL_M != QGEMM_DEFAULT_UNROLL_N
  .qtrmm_iunucopy = qtrmm_iunucopyTS,
  .qtrmm_iunncopy = qtrmm_iunncopyTS,
  .qtrmm_iutucopy = qtrmm_iutucopyTS,
  .qtrmm_iutncopy = qtrmm_iutncopyTS,
  .qtrmm_ilnucopy = qtrmm_ilnucopyTS,
  .qtrmm_ilnncopy = qtrmm_ilnncopyTS,
  .qtrmm_iltucopy = qtrmm_iltucopyTS,
  .qtrmm_iltncopy = qtrmm_iltncopyTS,
#else
  .qtrmm_iunucopy = qtrmm_ounucopyTS,
  .qtrmm_iunncopy = qtrmm_ounncopyTS,
  .qtrmm_iutucopy = qtrmm_outucopyTS,
  .qtrmm_iutncopy = qtrmm_outncopyTS,
  .qtrmm_ilnucopy = qtrmm_olnucopyTS,
  .qtrmm_ilnncopy = qtrmm_olnncopyTS,
  .qtrmm_iltucopy = qtrmm_oltucopyTS,
  .qtrmm_iltncopy = qtrmm_oltncopyTS,
#endif
  .qtrmm_ounucopy = qtrmm_ounucopyTS,
  .qtrmm_ounncopy = qtrmm_ounncopyTS,
  .qtrmm_outucopy = qtrmm_outucopyTS,
  .qtrmm_outncopy = qtrmm_outncopyTS,
  .qtrmm_olnucopy = qtrmm_olnucopyTS,
  .qtrmm_olnncopy = qtrmm_olnncopyTS,
  .qtrmm_oltucopy = qtrmm_oltucopyTS,
  .qtrmm_oltncopy = qtrmm_oltncopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_qsymm_dispatch_t openblas_qsymm_dispatchTS = {
#if QGEMM_DEFAULT_UNROLL_M != QGEMM_DEFAULT_UNROLL_N
  .qsymm_iutcopy = qsymm_iutcopyTS,
  .qsymm_iltcopy = qsymm_iltcopyTS,
#else
  .qsymm_iutcopy = qsymm_outcopyTS,
  .qsymm_iltcopy = qsymm_oltcopyTS,
#endif
  .qsymm_outcopy = qsymm_outcopyTS,
  .qsymm_oltcopy = qsymm_oltcopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_qneg_dispatch_t openblas_qneg_dispatchTS = {
#ifndef NO_LAPACK
  .qneg_tcopy = qneg_tcopyTS,
#else
  .qneg_tcopy = NULL,
#endif
};
#endif

#ifdef EXPRECISION
const openblas_qlaswp_dispatch_t openblas_qlaswp_dispatchTS = {
#ifndef NO_LAPACK
  .qlaswp_ncopy = qlaswp_ncopyTS,
#else
  .qlaswp_ncopy = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_camax_dispatch_t openblas_camax_dispatchTS = {
  .camax_k = camax_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_camin_dispatch_t openblas_camin_dispatchTS = {
  .camin_k = camin_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_icamax_dispatch_t openblas_icamax_dispatchTS = {
  .icamax_k = icamax_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_icamin_dispatch_t openblas_icamin_dispatchTS = {
  .icamin_k = icamin_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cnrm2_dispatch_t openblas_cnrm2_dispatchTS = {
  .cnrm2_k = cnrm2_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_casum_dispatch_t openblas_casum_dispatchTS = {
  .casum_k = casum_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_csum_dispatch_t openblas_csum_dispatchTS = {
  .csum_k = csum_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_ccopy_dispatch_t openblas_ccopy_dispatchTS = {
  .ccopy_k = ccopy_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cdotu_dispatch_t openblas_cdotu_dispatchTS = {
  .cdotu_k = cdotu_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cdotc_dispatch_t openblas_cdotc_dispatchTS = {
  .cdotc_k = cdotc_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_csrot_dispatch_t openblas_csrot_dispatchTS = {
 .csrot_k = csrot_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_caxpy_dispatch_t openblas_caxpy_dispatchTS = {
  .caxpy_k = caxpy_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_caxpyc_dispatch_t openblas_caxpyc_dispatchTS = {
  .caxpyc_k = caxpyc_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cscal_dispatch_t openblas_cscal_dispatchTS = {
  .cscal_k = cscal_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cswap_dispatch_t openblas_cswap_dispatchTS = {
  .cswap_k = cswap_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgemv_dispatch_t openblas_cgemv_dispatchTS = {
  .cgemv_n = cgemv_nTS,
  .cgemv_t = cgemv_tTS,
  .cgemv_r = cgemv_rTS,
  .cgemv_c = cgemv_cTS,
  .cgemv_o = cgemv_oTS,
  .cgemv_u = cgemv_uTS,
  .cgemv_s = cgemv_sTS,
  .cgemv_d = cgemv_dTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgeru_dispatch_t openblas_cgeru_dispatchTS = {
  .cgeru_k = cgeru_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgerc_dispatch_t openblas_cgerc_dispatchTS = {
  .cgerc_k = cgerc_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgerv_dispatch_t openblas_cgerv_dispatchTS = {
  .cgerv_k = cgerv_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgerd_dispatch_t openblas_cgerd_dispatchTS = {
  .cgerd_k = cgerd_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_csymv_dispatch_t openblas_csymv_dispatchTS = {
  .csymv_L = csymv_LTS,
  .csymv_U = csymv_UTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_chemv_dispatch_t openblas_chemv_dispatchTS = {
  .chemv_L = chemv_LTS,
  .chemv_U = chemv_UTS,
  .chemv_M = chemv_MTS,
  .chemv_V = chemv_VTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgemm_dispatch_t openblas_cgemm_dispatchTS = {
#ifdef ARCH_ARM64
#ifdef HAVE_SME
  .sme_cgemm_kernel = sme_cgemm_kernelTS,
#else
  .sme_cgemm_kernel = NULL,
#endif
#endif
  .cgemm_kernel_n = cgemm_kernel_nTS,
  .cgemm_kernel_l = cgemm_kernel_lTS,
  .cgemm_kernel_r = cgemm_kernel_rTS,
  .cgemm_kernel_b = cgemm_kernel_bTS,
  .cgemm_beta = cgemm_betaTS,
#if CGEMM_DEFAULT_UNROLL_M != CGEMM_DEFAULT_UNROLL_N
  .cgemm_incopy = cgemm_incopyTS,
  .cgemm_itcopy = cgemm_itcopyTS,
#else
  .cgemm_incopy = cgemm_oncopyTS,
  .cgemm_itcopy = cgemm_otcopyTS,
#endif
  .cgemm_oncopy = cgemm_oncopyTS,
  .cgemm_otcopy = cgemm_otcopyTS,
#ifdef SMALL_MATRIX_OPT
  .cgemm_small_matrix_permit = cgemm_small_matrix_permitTS,
  .cgemm_small_kernel_nn = cgemm_small_kernel_nnTS,
  .cgemm_small_kernel_nt = cgemm_small_kernel_ntTS,
  .cgemm_small_kernel_nr = cgemm_small_kernel_nrTS,
  .cgemm_small_kernel_nc = cgemm_small_kernel_ncTS,
  .cgemm_small_kernel_tn = cgemm_small_kernel_tnTS,
  .cgemm_small_kernel_tt = cgemm_small_kernel_ttTS,
  .cgemm_small_kernel_tr = cgemm_small_kernel_trTS,
  .cgemm_small_kernel_tc = cgemm_small_kernel_tcTS,
  .cgemm_small_kernel_rn = cgemm_small_kernel_rnTS,
  .cgemm_small_kernel_rt = cgemm_small_kernel_rtTS,
  .cgemm_small_kernel_rr = cgemm_small_kernel_rrTS,
  .cgemm_small_kernel_rc = cgemm_small_kernel_rcTS,
  .cgemm_small_kernel_cn = cgemm_small_kernel_cnTS,
  .cgemm_small_kernel_ct = cgemm_small_kernel_ctTS,
  .cgemm_small_kernel_cr = cgemm_small_kernel_crTS,
  .cgemm_small_kernel_cc = cgemm_small_kernel_ccTS,
  .cgemm_small_kernel_b0_nn = cgemm_small_kernel_b0_nnTS,
  .cgemm_small_kernel_b0_nt = cgemm_small_kernel_b0_ntTS,
  .cgemm_small_kernel_b0_nr = cgemm_small_kernel_b0_nrTS,
  .cgemm_small_kernel_b0_nc = cgemm_small_kernel_b0_ncTS,
  .cgemm_small_kernel_b0_tn = cgemm_small_kernel_b0_tnTS,
  .cgemm_small_kernel_b0_tt = cgemm_small_kernel_b0_ttTS,
  .cgemm_small_kernel_b0_tr = cgemm_small_kernel_b0_trTS,
  .cgemm_small_kernel_b0_tc = cgemm_small_kernel_b0_tcTS,
  .cgemm_small_kernel_b0_rn = cgemm_small_kernel_b0_rnTS,
  .cgemm_small_kernel_b0_rt = cgemm_small_kernel_b0_rtTS,
  .cgemm_small_kernel_b0_rr = cgemm_small_kernel_b0_rrTS,
  .cgemm_small_kernel_b0_rc = cgemm_small_kernel_b0_rcTS,
  .cgemm_small_kernel_b0_cn = cgemm_small_kernel_b0_cnTS,
  .cgemm_small_kernel_b0_ct = cgemm_small_kernel_b0_ctTS,
  .cgemm_small_kernel_b0_cr = cgemm_small_kernel_b0_crTS,
  .cgemm_small_kernel_b0_cc = cgemm_small_kernel_b0_ccTS,
#endif
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_csymm_dispatch_t openblas_csymm_dispatchTS = {
  .csymm_kernel_n = csymm_kernel_nTS,
  .csymm_kernel_l = csymm_kernel_lTS,
  .csymm_kernel_r = csymm_kernel_rTS,
  .csymm_kernel_b = csymm_kernel_bTS,
  .csymm_incopy = csymm_incopyTS,
  .csymm_itcopy = csymm_itcopyTS,
#if CGEMM_DEFAULT_UNROLL_M != CGEMM_DEFAULT_UNROLL_N
  .csymm_iutcopy = csymm_iutcopyTS,
  .csymm_iltcopy = csymm_iltcopyTS,
#else
  .csymm_iutcopy = csymm_outcopyTS,
  .csymm_iltcopy = csymm_oltcopyTS,
#endif
  .csymm_outcopy = csymm_outcopyTS,
  .csymm_oltcopy = csymm_oltcopyTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_ctrmm_dispatch_t openblas_ctrmm_dispatchTS = {
  .ctrmm_gemm_kernel_n = ctrmm_gemm_kernel_nTS,
  .ctrmm_gemm_kernel_l = ctrmm_gemm_kernel_lTS,
  .ctrmm_gemm_kernel_r = ctrmm_gemm_kernel_rTS,
  .ctrmm_gemm_kernel_b = ctrmm_gemm_kernel_bTS,
  .ctrmm_kernel_RN = ctrmm_kernel_RNTS,
  .ctrmm_kernel_RT = ctrmm_kernel_RTTS,
  .ctrmm_kernel_RR = ctrmm_kernel_RRTS,
  .ctrmm_kernel_RC = ctrmm_kernel_RCTS,
  .ctrmm_kernel_LN = ctrmm_kernel_LNTS,
  .ctrmm_kernel_LT = ctrmm_kernel_LTTS,
  .ctrmm_kernel_LR = ctrmm_kernel_LRTS,
  .ctrmm_kernel_LC = ctrmm_kernel_LCTS,
#if CGEMM_DEFAULT_UNROLL_M != CGEMM_DEFAULT_UNROLL_N
  .ctrmm_iunucopy = ctrmm_iunucopyTS,
  .ctrmm_iunncopy = ctrmm_iunncopyTS,
  .ctrmm_iutucopy = ctrmm_iutucopyTS,
  .ctrmm_iutncopy = ctrmm_iutncopyTS,
  .ctrmm_ilnucopy = ctrmm_ilnucopyTS,
  .ctrmm_ilnncopy = ctrmm_ilnncopyTS,
  .ctrmm_iltucopy = ctrmm_iltucopyTS,
  .ctrmm_iltncopy = ctrmm_iltncopyTS,
#else
  .ctrmm_iunucopy = ctrmm_ounucopyTS,
  .ctrmm_iunncopy = ctrmm_ounncopyTS,
  .ctrmm_iutucopy = ctrmm_outucopyTS,
  .ctrmm_iutncopy = ctrmm_outncopyTS,
  .ctrmm_ilnucopy = ctrmm_olnucopyTS,
  .ctrmm_ilnncopy = ctrmm_olnncopyTS,
  .ctrmm_iltucopy = ctrmm_oltucopyTS,
  .ctrmm_iltncopy = ctrmm_oltncopyTS,
#endif
  .ctrmm_incopy = ctrmm_incopyTS,
  .ctrmm_itcopy = ctrmm_itcopyTS,
  .ctrmm_ounucopy = ctrmm_ounucopyTS,
  .ctrmm_ounncopy = ctrmm_ounncopyTS,
  .ctrmm_outucopy = ctrmm_outucopyTS,
  .ctrmm_outncopy = ctrmm_outncopyTS,
  .ctrmm_olnucopy = ctrmm_olnucopyTS,
  .ctrmm_olnncopy = ctrmm_olnncopyTS,
  .ctrmm_oltucopy = ctrmm_oltucopyTS,
  .ctrmm_oltncopy = ctrmm_oltncopyTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_ctrsm_dispatch_t openblas_ctrsm_dispatchTS = {
  .ctrsm_kernel_LN = ctrsm_kernel_LNTS,
  .ctrsm_kernel_LT = ctrsm_kernel_LTTS,
  .ctrsm_kernel_LR = ctrsm_kernel_LRTS,
  .ctrsm_kernel_LC = ctrsm_kernel_LCTS,
  .ctrsm_kernel_RN = ctrsm_kernel_RNTS,
  .ctrsm_kernel_RT = ctrsm_kernel_RTTS,
  .ctrsm_kernel_RR = ctrsm_kernel_RRTS,
  .ctrsm_kernel_RC = ctrsm_kernel_RCTS,
#if CGEMM_DEFAULT_UNROLL_M != CGEMM_DEFAULT_UNROLL_N
  .ctrsm_iunucopy = ctrsm_iunucopyTS,
  .ctrsm_iunncopy = ctrsm_iunncopyTS,
  .ctrsm_iutucopy = ctrsm_iutucopyTS,
  .ctrsm_iutncopy = ctrsm_iutncopyTS,
  .ctrsm_ilnucopy = ctrsm_ilnucopyTS,
  .ctrsm_ilnncopy = ctrsm_ilnncopyTS,
  .ctrsm_iltucopy = ctrsm_iltucopyTS,
  .ctrsm_iltncopy = ctrsm_iltncopyTS,
#else
  .ctrsm_iunucopy = ctrsm_ounucopyTS,
  .ctrsm_iunncopy = ctrsm_ounncopyTS,
  .ctrsm_iutucopy = ctrsm_outucopyTS,
  .ctrsm_iutncopy = ctrsm_outncopyTS,
  .ctrsm_ilnucopy = ctrsm_olnucopyTS,
  .ctrsm_ilnncopy = ctrsm_olnncopyTS,
  .ctrsm_iltucopy = ctrsm_oltucopyTS,
  .ctrsm_iltncopy = ctrsm_oltncopyTS,
#endif
  .ctrsm_ounucopy = ctrsm_ounucopyTS,
  .ctrsm_ounncopy = ctrsm_ounncopyTS,
  .ctrsm_outucopy = ctrsm_outucopyTS,
  .ctrsm_outncopy = ctrsm_outncopyTS,
  .ctrsm_olnucopy = ctrsm_olnucopyTS,
  .ctrsm_olnncopy = ctrsm_olnncopyTS,
  .ctrsm_oltucopy = ctrsm_oltucopyTS,
  .ctrsm_oltncopy = ctrsm_oltncopyTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_chemm_dispatch_t openblas_chemm_dispatchTS = {
#if CGEMM_DEFAULT_UNROLL_M != CGEMM_DEFAULT_UNROLL_N
  .chemm_iutcopy = chemm_iutcopyTS,
  .chemm_iltcopy = chemm_iltcopyTS,
#else
  .chemm_iutcopy = chemm_outcopyTS,
  .chemm_iltcopy = chemm_oltcopyTS,
#endif
  .chemm_outcopy = chemm_outcopyTS,
  .chemm_oltcopy = chemm_oltcopyTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgemm3m_dispatch_t openblas_cgemm3m_dispatchTS = {
#if (USE_GEMM3M)
  .cgemm3m_kernel = cgemm3m_kernelTS,
  .cgemm3m_incopyb = cgemm3m_incopybTS,
  .cgemm3m_incopyr = cgemm3m_incopyrTS,
  .cgemm3m_incopyi = cgemm3m_incopyiTS,
  .cgemm3m_itcopyb = cgemm3m_itcopybTS,
  .cgemm3m_itcopyr = cgemm3m_itcopyrTS,
  .cgemm3m_itcopyi = cgemm3m_itcopyiTS,
  .cgemm3m_oncopyb = cgemm3m_oncopybTS,
  .cgemm3m_oncopyr = cgemm3m_oncopyrTS,
  .cgemm3m_oncopyi = cgemm3m_oncopyiTS,
  .cgemm3m_otcopyb = cgemm3m_otcopybTS,
  .cgemm3m_otcopyr = cgemm3m_otcopyrTS,
  .cgemm3m_otcopyi = cgemm3m_otcopyiTS,
#else
  .cgemm3m_kernel = NULL,
  .cgemm3m_incopyb = NULL,
  .cgemm3m_incopyr = NULL,
  .cgemm3m_incopyi = NULL,
  .cgemm3m_itcopyb = NULL,
  .cgemm3m_itcopyr = NULL,
  .cgemm3m_itcopyi = NULL,
  .cgemm3m_oncopyb = NULL,
  .cgemm3m_oncopyr = NULL,
  .cgemm3m_oncopyi = NULL,
  .cgemm3m_otcopyb = NULL,
  .cgemm3m_otcopyr = NULL,
  .cgemm3m_otcopyi = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_csymm3m_dispatch_t openblas_csymm3m_dispatchTS = {
#if (USE_GEMM3M)
  .csymm3m_iucopyb = csymm3m_iucopybTS,
  .csymm3m_ilcopyb = csymm3m_ilcopybTS,
  .csymm3m_iucopyr = csymm3m_iucopyrTS,
  .csymm3m_ilcopyr = csymm3m_ilcopyrTS,
  .csymm3m_iucopyi = csymm3m_iucopyiTS,
  .csymm3m_ilcopyi = csymm3m_ilcopyiTS,
  .csymm3m_oucopyb = csymm3m_oucopybTS,
  .csymm3m_olcopyb = csymm3m_olcopybTS,
  .csymm3m_oucopyr = csymm3m_oucopyrTS,
  .csymm3m_olcopyr = csymm3m_olcopyrTS,
  .csymm3m_oucopyi = csymm3m_oucopyiTS,
  .csymm3m_olcopyi = csymm3m_olcopyiTS,
#else
  .csymm3m_iucopyb = NULL,
  .csymm3m_ilcopyb = NULL,
  .csymm3m_iucopyr = NULL,
  .csymm3m_ilcopyr = NULL,
  .csymm3m_iucopyi = NULL,
  .csymm3m_ilcopyi = NULL,
  .csymm3m_oucopyb = NULL,
  .csymm3m_olcopyb = NULL,
  .csymm3m_oucopyr = NULL,
  .csymm3m_olcopyr = NULL,
  .csymm3m_oucopyi = NULL,
  .csymm3m_olcopyi = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_chemm3m_dispatch_t openblas_chemm3m_dispatchTS = {
#if (USE_GEMM3M)
  .chemm3m_iucopyb = chemm3m_iucopybTS,
  .chemm3m_ilcopyb = chemm3m_ilcopybTS,
  .chemm3m_iucopyr = chemm3m_iucopyrTS,
  .chemm3m_ilcopyr = chemm3m_ilcopyrTS,
  .chemm3m_iucopyi = chemm3m_iucopyiTS,
  .chemm3m_ilcopyi = chemm3m_ilcopyiTS,
  .chemm3m_oucopyb = chemm3m_oucopybTS,
  .chemm3m_olcopyb = chemm3m_olcopybTS,
  .chemm3m_oucopyr = chemm3m_oucopyrTS,
  .chemm3m_olcopyr = chemm3m_olcopyrTS,
  .chemm3m_oucopyi = chemm3m_oucopyiTS,
  .chemm3m_olcopyi = chemm3m_olcopyiTS,
#else
  .chemm3m_iucopyb = NULL,
  .chemm3m_ilcopyb = NULL,
  .chemm3m_iucopyr = NULL,
  .chemm3m_ilcopyr = NULL,
  .chemm3m_iucopyi = NULL,
  .chemm3m_ilcopyi = NULL,
  .chemm3m_oucopyb = NULL,
  .chemm3m_olcopyb = NULL,
  .chemm3m_oucopyr = NULL,
  .chemm3m_olcopyr = NULL,
  .chemm3m_oucopyi = NULL,
  .chemm3m_olcopyi = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cneg_dispatch_t openblas_cneg_dispatchTS = {
#ifndef NO_LAPACK
  .cneg_tcopy = cneg_tcopyTS,
#else
  .cneg_tcopy = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_claswp_dispatch_t openblas_claswp_dispatchTS = {
#ifndef NO_LAPACK
   .claswp_ncopy = claswp_ncopyTS,
#else
  .claswp_ncopy = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zamax_dispatch_t openblas_zamax_dispatchTS = {
  .zamax_k = zamax_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zamin_dispatch_t openblas_zamin_dispatchTS = {
  .zamin_k = zamin_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_izamax_dispatch_t openblas_izamax_dispatchTS = {
  .izamax_k = izamax_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_izamin_dispatch_t openblas_izamin_dispatchTS = {
  .izamin_k = izamin_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_znrm2_dispatch_t openblas_znrm2_dispatchTS = {
  .znrm2_k = znrm2_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zasum_dispatch_t openblas_zasum_dispatchTS = {
  .zasum_k = zasum_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zsum_dispatch_t openblas_zsum_dispatchTS = {
  .zsum_k = zsum_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zcopy_dispatch_t openblas_zcopy_dispatchTS = {
  .zcopy_k = zcopy_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zdotu_dispatch_t openblas_zdotu_dispatchTS = {
  .zdotu_k = zdotu_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zdotc_dispatch_t openblas_zdotc_dispatchTS = {
  .zdotc_k = zdotc_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zdrot_dispatch_t openblas_zdrot_dispatchTS = {
  .zdrot_k = zdrot_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zaxpy_dispatch_t openblas_zaxpy_dispatchTS = {
  .zaxpy_k = zaxpy_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zaxpyc_dispatch_t openblas_zaxpyc_dispatchTS = {
  .zaxpyc_k = zaxpyc_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zscal_dispatch_t openblas_zscal_dispatchTS = {
  .zscal_k = zscal_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zswap_dispatch_t openblas_zswap_dispatchTS = {
  .zswap_k = zswap_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zgemv_dispatch_t openblas_zgemv_dispatchTS = {
  .zgemv_n = zgemv_nTS,
  .zgemv_t = zgemv_tTS,
  .zgemv_r = zgemv_rTS,
  .zgemv_c = zgemv_cTS,
  .zgemv_o = zgemv_oTS,
  .zgemv_u = zgemv_uTS,
  .zgemv_s = zgemv_sTS,
  .zgemv_d = zgemv_dTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zgeru_dispatch_t openblas_zgeru_dispatchTS = {
  .zgeru_k = zgeru_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zgerc_dispatch_t openblas_zgerc_dispatchTS = {
  .zgerc_k = zgerc_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zgerv_dispatch_t openblas_zgerv_dispatchTS = {
  .zgerv_k = zgerv_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zgerd_dispatch_t openblas_zgerd_dispatchTS = {
  .zgerd_k = zgerd_kTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zsymv_dispatch_t openblas_zsymv_dispatchTS = {
  .zsymv_L = zsymv_LTS,
  .zsymv_U = zsymv_UTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zhemv_dispatch_t openblas_zhemv_dispatchTS = {
  .zhemv_L = zhemv_LTS,
  .zhemv_U = zhemv_UTS,
  .zhemv_M = zhemv_MTS,
  .zhemv_V = zhemv_VTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zgemm_dispatch_t openblas_zgemm_dispatchTS = {
#ifdef ARCH_ARM64
#ifdef HAVE_SME
  .sme_zgemm_kernel = sme_zgemm_kernelTS,
#else
  .sme_zgemm_kernel = NULL,
#endif
#endif
  .zgemm_kernel_n = zgemm_kernel_nTS,
  .zgemm_kernel_l = zgemm_kernel_lTS,
  .zgemm_kernel_r = zgemm_kernel_rTS,
  .zgemm_kernel_b = zgemm_kernel_bTS,
  .zgemm_beta = zgemm_betaTS,
#if ZGEMM_DEFAULT_UNROLL_M != ZGEMM_DEFAULT_UNROLL_N
  .zgemm_incopy = zgemm_incopyTS,
  .zgemm_itcopy = zgemm_itcopyTS,
#else
  .zgemm_incopy = zgemm_oncopyTS,
  .zgemm_itcopy = zgemm_otcopyTS,
#endif
  .zgemm_oncopy = zgemm_oncopyTS,
  .zgemm_otcopy = zgemm_otcopyTS,
#ifdef SMALL_MATRIX_OPT
  .zgemm_small_matrix_permit = zgemm_small_matrix_permitTS,
  .zgemm_small_kernel_nn = zgemm_small_kernel_nnTS,
  .zgemm_small_kernel_nt = zgemm_small_kernel_ntTS,
  .zgemm_small_kernel_nr = zgemm_small_kernel_nrTS,
  .zgemm_small_kernel_nc = zgemm_small_kernel_ncTS,
  .zgemm_small_kernel_tn = zgemm_small_kernel_tnTS,
  .zgemm_small_kernel_tt = zgemm_small_kernel_ttTS,
  .zgemm_small_kernel_tr = zgemm_small_kernel_trTS,
  .zgemm_small_kernel_tc = zgemm_small_kernel_tcTS,
  .zgemm_small_kernel_rn = zgemm_small_kernel_rnTS,
  .zgemm_small_kernel_rt = zgemm_small_kernel_rtTS,
  .zgemm_small_kernel_rr = zgemm_small_kernel_rrTS,
  .zgemm_small_kernel_rc = zgemm_small_kernel_rcTS,
  .zgemm_small_kernel_cn = zgemm_small_kernel_cnTS,
  .zgemm_small_kernel_ct = zgemm_small_kernel_ctTS,
  .zgemm_small_kernel_cr = zgemm_small_kernel_crTS,
  .zgemm_small_kernel_cc = zgemm_small_kernel_ccTS,
  .zgemm_small_kernel_b0_nn = zgemm_small_kernel_b0_nnTS,
  .zgemm_small_kernel_b0_nt = zgemm_small_kernel_b0_ntTS,
  .zgemm_small_kernel_b0_nr = zgemm_small_kernel_b0_nrTS,
  .zgemm_small_kernel_b0_nc = zgemm_small_kernel_b0_ncTS,
  .zgemm_small_kernel_b0_tn = zgemm_small_kernel_b0_tnTS,
  .zgemm_small_kernel_b0_tt = zgemm_small_kernel_b0_ttTS,
  .zgemm_small_kernel_b0_tr = zgemm_small_kernel_b0_trTS,
  .zgemm_small_kernel_b0_tc = zgemm_small_kernel_b0_tcTS,
  .zgemm_small_kernel_b0_rn = zgemm_small_kernel_b0_rnTS,
  .zgemm_small_kernel_b0_rt = zgemm_small_kernel_b0_rtTS,
  .zgemm_small_kernel_b0_rr = zgemm_small_kernel_b0_rrTS,
  .zgemm_small_kernel_b0_rc = zgemm_small_kernel_b0_rcTS,
  .zgemm_small_kernel_b0_cn = zgemm_small_kernel_b0_cnTS,
  .zgemm_small_kernel_b0_ct = zgemm_small_kernel_b0_ctTS,
  .zgemm_small_kernel_b0_cr = zgemm_small_kernel_b0_crTS,
  .zgemm_small_kernel_b0_cc = zgemm_small_kernel_b0_ccTS,
#endif
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zsymm_dispatch_t openblas_zsymm_dispatchTS = {
  .zsymm_kernel_n = zsymm_kernel_nTS,
  .zsymm_kernel_l = zsymm_kernel_lTS,
  .zsymm_kernel_r = zsymm_kernel_rTS,
  .zsymm_kernel_b = zsymm_kernel_bTS,
  .zsymm_incopy = zsymm_incopyTS,
  .zsymm_itcopy = zsymm_itcopyTS,
#if ZGEMM_DEFAULT_UNROLL_M != ZGEMM_DEFAULT_UNROLL_N
  .zsymm_iutcopy = zsymm_iutcopyTS,
  .zsymm_iltcopy = zsymm_iltcopyTS,
#else
  .zsymm_iutcopy = zsymm_outcopyTS,
  .zsymm_iltcopy = zsymm_oltcopyTS,
#endif
  .zsymm_outcopy = zsymm_outcopyTS,
  .zsymm_oltcopy = zsymm_oltcopyTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_ztrmm_dispatch_t openblas_ztrmm_dispatchTS = {
  .ztrmm_gemm_kernel_n = ztrmm_gemm_kernel_nTS,
  .ztrmm_gemm_kernel_l = ztrmm_gemm_kernel_lTS,
  .ztrmm_gemm_kernel_r = ztrmm_gemm_kernel_rTS,
  .ztrmm_gemm_kernel_b = ztrmm_gemm_kernel_bTS,
  .ztrmm_kernel_RN = ztrmm_kernel_RNTS,
  .ztrmm_kernel_RT = ztrmm_kernel_RTTS,
  .ztrmm_kernel_RR = ztrmm_kernel_RRTS,
  .ztrmm_kernel_RC = ztrmm_kernel_RCTS,
  .ztrmm_kernel_LN = ztrmm_kernel_LNTS,
  .ztrmm_kernel_LT = ztrmm_kernel_LTTS,
  .ztrmm_kernel_LR = ztrmm_kernel_LRTS,
  .ztrmm_kernel_LC = ztrmm_kernel_LCTS,
#if ZGEMM_DEFAULT_UNROLL_M != ZGEMM_DEFAULT_UNROLL_N
  .ztrmm_iunucopy = ztrmm_iunucopyTS,
  .ztrmm_iunncopy = ztrmm_iunncopyTS,
  .ztrmm_iutucopy = ztrmm_iutucopyTS,
  .ztrmm_iutncopy = ztrmm_iutncopyTS,
  .ztrmm_ilnucopy = ztrmm_ilnucopyTS,
  .ztrmm_ilnncopy = ztrmm_ilnncopyTS,
  .ztrmm_iltucopy = ztrmm_iltucopyTS,
  .ztrmm_iltncopy = ztrmm_iltncopyTS,
#else
  .ztrmm_iunucopy = ztrmm_ounucopyTS,
  .ztrmm_iunncopy = ztrmm_ounncopyTS,
  .ztrmm_iutucopy = ztrmm_outucopyTS,
  .ztrmm_iutncopy = ztrmm_outncopyTS,
  .ztrmm_ilnucopy = ztrmm_olnucopyTS,
  .ztrmm_ilnncopy = ztrmm_olnncopyTS,
  .ztrmm_iltucopy = ztrmm_oltucopyTS,
  .ztrmm_iltncopy = ztrmm_oltncopyTS,
#endif
  .ztrmm_incopy = ztrmm_incopyTS,
  .ztrmm_itcopy = ztrmm_itcopyTS,
  .ztrmm_ounucopy = ztrmm_ounucopyTS,
  .ztrmm_ounncopy = ztrmm_ounncopyTS,
  .ztrmm_outucopy = ztrmm_outucopyTS,
  .ztrmm_outncopy = ztrmm_outncopyTS,
  .ztrmm_olnucopy = ztrmm_olnucopyTS,
  .ztrmm_olnncopy = ztrmm_olnncopyTS,
  .ztrmm_oltucopy = ztrmm_oltucopyTS,
  .ztrmm_oltncopy = ztrmm_oltncopyTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_ztrsm_dispatch_t openblas_ztrsm_dispatchTS = {
  .ztrsm_kernel_LN = ztrsm_kernel_LNTS,
  .ztrsm_kernel_LT = ztrsm_kernel_LTTS,
  .ztrsm_kernel_LR = ztrsm_kernel_LRTS,
  .ztrsm_kernel_LC = ztrsm_kernel_LCTS,
  .ztrsm_kernel_RN = ztrsm_kernel_RNTS,
  .ztrsm_kernel_RT = ztrsm_kernel_RTTS,
  .ztrsm_kernel_RR = ztrsm_kernel_RRTS,
  .ztrsm_kernel_RC = ztrsm_kernel_RCTS,
#if ZGEMM_DEFAULT_UNROLL_M != ZGEMM_DEFAULT_UNROLL_N
  .ztrsm_iunucopy = ztrsm_iunucopyTS,
  .ztrsm_iunncopy = ztrsm_iunncopyTS,
  .ztrsm_iutucopy = ztrsm_iutucopyTS,
  .ztrsm_iutncopy = ztrsm_iutncopyTS,
  .ztrsm_ilnucopy = ztrsm_ilnucopyTS,
  .ztrsm_ilnncopy = ztrsm_ilnncopyTS,
  .ztrsm_iltucopy = ztrsm_iltucopyTS,
  .ztrsm_iltncopy = ztrsm_iltncopyTS,
#else
  .ztrsm_iunucopy = ztrsm_ounucopyTS,
  .ztrsm_iunncopy = ztrsm_ounncopyTS,
  .ztrsm_iutucopy = ztrsm_outucopyTS,
  .ztrsm_iutncopy = ztrsm_outncopyTS,
  .ztrsm_ilnucopy = ztrsm_olnucopyTS,
  .ztrsm_ilnncopy = ztrsm_olnncopyTS,
  .ztrsm_iltucopy = ztrsm_oltucopyTS,
  .ztrsm_iltncopy = ztrsm_oltncopyTS,
#endif
  .ztrsm_ounucopy = ztrsm_ounucopyTS,
  .ztrsm_ounncopy = ztrsm_ounncopyTS,
  .ztrsm_outucopy = ztrsm_outucopyTS,
  .ztrsm_outncopy = ztrsm_outncopyTS,
  .ztrsm_olnucopy = ztrsm_olnucopyTS,
  .ztrsm_olnncopy = ztrsm_olnncopyTS,
  .ztrsm_oltucopy = ztrsm_oltucopyTS,
  .ztrsm_oltncopy = ztrsm_oltncopyTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zhemm_dispatch_t openblas_zhemm_dispatchTS = {
#if ZGEMM_DEFAULT_UNROLL_M != ZGEMM_DEFAULT_UNROLL_N
  .zhemm_iutcopy = zhemm_iutcopyTS,
  .zhemm_iltcopy = zhemm_iltcopyTS,
#else
  .zhemm_iutcopy = zhemm_outcopyTS,
  .zhemm_iltcopy = zhemm_oltcopyTS,
#endif
  .zhemm_outcopy = zhemm_outcopyTS,
  .zhemm_oltcopy = zhemm_oltcopyTS,
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zgemm3m_dispatch_t openblas_zgemm3m_dispatchTS = {
#if (USE_GEMM3M)
  .zgemm3m_kernel = zgemm3m_kernelTS,
  .zgemm3m_incopyb = zgemm3m_incopybTS,
  .zgemm3m_incopyr = zgemm3m_incopyrTS,
  .zgemm3m_incopyi = zgemm3m_incopyiTS,
  .zgemm3m_itcopyb = zgemm3m_itcopybTS,
  .zgemm3m_itcopyr = zgemm3m_itcopyrTS,
  .zgemm3m_itcopyi = zgemm3m_itcopyiTS,
  .zgemm3m_oncopyb = zgemm3m_oncopybTS,
  .zgemm3m_oncopyr = zgemm3m_oncopyrTS,
  .zgemm3m_oncopyi = zgemm3m_oncopyiTS,
  .zgemm3m_otcopyb = zgemm3m_otcopybTS,
  .zgemm3m_otcopyr = zgemm3m_otcopyrTS,
  .zgemm3m_otcopyi = zgemm3m_otcopyiTS,
#else
  .zgemm3m_kernel = NULL,
  .zgemm3m_incopyb = NULL,
  .zgemm3m_incopyr = NULL,
  .zgemm3m_incopyi = NULL,
  .zgemm3m_itcopyb = NULL,
  .zgemm3m_itcopyr = NULL,
  .zgemm3m_itcopyi = NULL,
  .zgemm3m_oncopyb = NULL,
  .zgemm3m_oncopyr = NULL,
  .zgemm3m_oncopyi = NULL,
  .zgemm3m_otcopyb = NULL,
  .zgemm3m_otcopyr = NULL,
  .zgemm3m_otcopyi = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zsymm3m_dispatch_t openblas_zsymm3m_dispatchTS = {
#if (USE_GEMM3M)
  .zsymm3m_iucopyb = zsymm3m_iucopybTS,
  .zsymm3m_ilcopyb = zsymm3m_ilcopybTS,
  .zsymm3m_iucopyr = zsymm3m_iucopyrTS,
  .zsymm3m_ilcopyr = zsymm3m_ilcopyrTS,
  .zsymm3m_iucopyi = zsymm3m_iucopyiTS,
  .zsymm3m_ilcopyi = zsymm3m_ilcopyiTS,
  .zsymm3m_oucopyb = zsymm3m_oucopybTS,
  .zsymm3m_olcopyb = zsymm3m_olcopybTS,
  .zsymm3m_oucopyr = zsymm3m_oucopyrTS,
  .zsymm3m_olcopyr = zsymm3m_olcopyrTS,
  .zsymm3m_oucopyi = zsymm3m_oucopyiTS,
  .zsymm3m_olcopyi = zsymm3m_olcopyiTS,
#else
  .zsymm3m_iucopyb = NULL,
  .zsymm3m_ilcopyb = NULL,
  .zsymm3m_iucopyr = NULL,
  .zsymm3m_ilcopyr = NULL,
  .zsymm3m_iucopyi = NULL,
  .zsymm3m_ilcopyi = NULL,
  .zsymm3m_oucopyb = NULL,
  .zsymm3m_olcopyb = NULL,
  .zsymm3m_oucopyr = NULL,
  .zsymm3m_olcopyr = NULL,
  .zsymm3m_oucopyi = NULL,
  .zsymm3m_olcopyi = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zhemm3m_dispatch_t openblas_zhemm3m_dispatchTS = {
#if (USE_GEMM3M)
  .zhemm3m_iucopyb = zhemm3m_iucopybTS,
  .zhemm3m_ilcopyb = zhemm3m_ilcopybTS,
  .zhemm3m_iucopyr = zhemm3m_iucopyrTS,
  .zhemm3m_ilcopyr = zhemm3m_ilcopyrTS,
  .zhemm3m_iucopyi = zhemm3m_iucopyiTS,
  .zhemm3m_ilcopyi = zhemm3m_ilcopyiTS,
  .zhemm3m_oucopyb = zhemm3m_oucopybTS,
  .zhemm3m_olcopyb = zhemm3m_olcopybTS,
  .zhemm3m_oucopyr = zhemm3m_oucopyrTS,
  .zhemm3m_olcopyr = zhemm3m_olcopyrTS,
  .zhemm3m_oucopyi = zhemm3m_oucopyiTS,
  .zhemm3m_olcopyi = zhemm3m_olcopyiTS,
#else
  .zhemm3m_iucopyb = NULL,
  .zhemm3m_ilcopyb = NULL,
  .zhemm3m_iucopyr = NULL,
  .zhemm3m_ilcopyr = NULL,
  .zhemm3m_iucopyi = NULL,
  .zhemm3m_ilcopyi = NULL,
  .zhemm3m_oucopyb = NULL,
  .zhemm3m_olcopyb = NULL,
  .zhemm3m_oucopyr = NULL,
  .zhemm3m_olcopyr = NULL,
  .zhemm3m_oucopyi = NULL,
  .zhemm3m_olcopyi = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zneg_dispatch_t openblas_zneg_dispatchTS = {
#ifndef NO_LAPACK
  .zneg_tcopy = zneg_tcopyTS,
#else
  .zneg_tcopy = NULL,
#endif
};
#endif

#if (BUILD_COMPLEX16 == 1)
const openblas_zlaswp_dispatch_t openblas_zlaswp_dispatchTS = {
#ifndef NO_LAPACK
  .zlaswp_ncopy = zlaswp_ncopyTS,
#else
  .zlaswp_ncopy = NULL,
#endif
};
#endif

#ifdef EXPRECISION
const openblas_xamax_dispatch_t openblas_xamax_dispatchTS = {
  .xamax_k = xamax_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xamin_dispatch_t openblas_xamin_dispatchTS = {
  .xamin_k = xamin_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_ixamax_dispatch_t openblas_ixamax_dispatchTS = {
  .ixamax_k = ixamax_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_ixamin_dispatch_t openblas_ixamin_dispatchTS = {
  .ixamin_k = ixamin_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xnrm2_dispatch_t openblas_xnrm2_dispatchTS = {
  .xnrm2_k = xnrm2_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xasum_dispatch_t openblas_xasum_dispatchTS = {
  .xasum_k = xasum_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xsum_dispatch_t openblas_xsum_dispatchTS = {
  .xsum_k = xsum_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xcopy_dispatch_t openblas_xcopy_dispatchTS = {
  .xcopy_k = xcopy_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xdotu_dispatch_t openblas_xdotu_dispatchTS = {
  .xdotu_k = xdotu_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xdotc_dispatch_t openblas_xdotc_dispatchTS = {
  .xdotc_k = xdotc_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xqrot_dispatch_t openblas_xqrot_dispatchTS = {
  .xqrot_k = xqrot_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xaxpy_dispatch_t openblas_xaxpy_dispatchTS = {
  .xaxpy_k = xaxpy_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xaxpyc_dispatch_t openblas_xaxpyc_dispatchTS = {
  .xaxpyc_k = xaxpyc_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xscal_dispatch_t openblas_xscal_dispatchTS = {
  .xscal_k = xscal_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xswap_dispatch_t openblas_xswap_dispatchTS = {
  .xswap_k = xswap_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xgemv_dispatch_t openblas_xgemv_dispatchTS = {
  .xgemv_n = xgemv_nTS,
  .xgemv_t = xgemv_tTS,
  .xgemv_r = xgemv_rTS,
  .xgemv_c = xgemv_cTS,
  .xgemv_o = xgemv_oTS,
  .xgemv_u = xgemv_uTS,
  .xgemv_s = xgemv_sTS,
  .xgemv_d = xgemv_dTS,
};
#endif

#ifdef EXPRECISION
const openblas_xgeru_dispatch_t openblas_xgeru_dispatchTS = {
  .xgeru_k = xgeru_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xgerc_dispatch_t openblas_xgerc_dispatchTS = {
  .xgerc_k = xgerc_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xgerv_dispatch_t openblas_xgerv_dispatchTS = {
  .xgerv_k = xgerv_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xgerd_dispatch_t openblas_xgerd_dispatchTS = {
  .xgerd_k = xgerd_kTS,
};
#endif

#ifdef EXPRECISION
const openblas_xsymv_dispatch_t openblas_xsymv_dispatchTS = {
  .xsymv_L = xsymv_LTS,
  .xsymv_U = xsymv_UTS,
};
#endif

#ifdef EXPRECISION
const openblas_xhemv_dispatch_t openblas_xhemv_dispatchTS = {
  .xhemv_L = xhemv_LTS,
  .xhemv_U = xhemv_UTS,
  .xhemv_M = xhemv_MTS,
  .xhemv_V = xhemv_VTS,
};
#endif

#ifdef EXPRECISION
const openblas_xgemm_dispatch_t openblas_xgemm_dispatchTS = {
  .xgemm_kernel_n = xgemm_kernel_nTS,
  .xgemm_kernel_l = xgemm_kernel_lTS,
  .xgemm_kernel_r = xgemm_kernel_rTS,
  .xgemm_kernel_b = xgemm_kernel_bTS,
  .xgemm_beta = xgemm_betaTS,
#if XGEMM_DEFAULT_UNROLL_M != XGEMM_DEFAULT_UNROLL_N
  .xgemm_incopy = xgemm_incopyTS,
  .xgemm_itcopy = xgemm_itcopyTS,
#else
  .xgemm_incopy = xgemm_oncopyTS,
  .xgemm_itcopy = xgemm_otcopyTS,
#endif
  .xgemm_oncopy = xgemm_oncopyTS,
  .xgemm_otcopy = xgemm_otcopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_xtrsm_dispatch_t openblas_xtrsm_dispatchTS = {
  .xtrsm_kernel_LN = xtrsm_kernel_LNTS,
  .xtrsm_kernel_LT = xtrsm_kernel_LTTS,
  .xtrsm_kernel_LR = xtrsm_kernel_LRTS,
  .xtrsm_kernel_LC = xtrsm_kernel_LCTS,
  .xtrsm_kernel_RN = xtrsm_kernel_RNTS,
  .xtrsm_kernel_RT = xtrsm_kernel_RTTS,
  .xtrsm_kernel_RR = xtrsm_kernel_RRTS,
  .xtrsm_kernel_RC = xtrsm_kernel_RCTS,
#if XGEMM_DEFAULT_UNROLL_M != XGEMM_DEFAULT_UNROLL_N
  .xtrsm_iunucopy = xtrsm_iunucopyTS,
  .xtrsm_iunncopy = xtrsm_iunncopyTS,
  .xtrsm_iutucopy = xtrsm_iutucopyTS,
  .xtrsm_iutncopy = xtrsm_iutncopyTS,
  .xtrsm_ilnucopy = xtrsm_ilnucopyTS,
  .xtrsm_ilnncopy = xtrsm_ilnncopyTS,
  .xtrsm_iltucopy = xtrsm_iltucopyTS,
  .xtrsm_iltncopy = xtrsm_iltncopyTS,
#else
  .xtrsm_iunucopy = xtrsm_ounucopyTS,
  .xtrsm_iunncopy = xtrsm_ounncopyTS,
  .xtrsm_iutucopy = xtrsm_outucopyTS,
  .xtrsm_iutncopy = xtrsm_outncopyTS,
  .xtrsm_ilnucopy = xtrsm_olnucopyTS,
  .xtrsm_ilnncopy = xtrsm_olnncopyTS,
  .xtrsm_iltucopy = xtrsm_oltucopyTS,
  .xtrsm_iltncopy = xtrsm_oltncopyTS,
#endif
  .xtrsm_ounucopy = xtrsm_ounucopyTS,
  .xtrsm_ounncopy = xtrsm_ounncopyTS,
  .xtrsm_outucopy = xtrsm_outucopyTS,
  .xtrsm_outncopy = xtrsm_outncopyTS,
  .xtrsm_olnucopy = xtrsm_olnucopyTS,
  .xtrsm_olnncopy = xtrsm_olnncopyTS,
  .xtrsm_oltucopy = xtrsm_oltucopyTS,
  .xtrsm_oltncopy = xtrsm_oltncopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_xtrmm_dispatch_t openblas_xtrmm_dispatchTS = {
  .xtrmm_kernel_RN = xtrmm_kernel_RNTS,
  .xtrmm_kernel_RT = xtrmm_kernel_RTTS,
  .xtrmm_kernel_RR = xtrmm_kernel_RRTS,
  .xtrmm_kernel_RC = xtrmm_kernel_RCTS,
  .xtrmm_kernel_LN = xtrmm_kernel_LNTS,
  .xtrmm_kernel_LT = xtrmm_kernel_LTTS,
  .xtrmm_kernel_LR = xtrmm_kernel_LRTS,
  .xtrmm_kernel_LC = xtrmm_kernel_LCTS,
#if XGEMM_DEFAULT_UNROLL_M != XGEMM_DEFAULT_UNROLL_N
  .xtrmm_iunucopy = xtrmm_iunucopyTS,
  .xtrmm_iunncopy = xtrmm_iunncopyTS,
  .xtrmm_iutucopy = xtrmm_iutucopyTS,
  .xtrmm_iutncopy = xtrmm_iutncopyTS,
  .xtrmm_ilnucopy = xtrmm_ilnucopyTS,
  .xtrmm_ilnncopy = xtrmm_ilnncopyTS,
  .xtrmm_iltucopy = xtrmm_iltucopyTS,
  .xtrmm_iltncopy = xtrmm_iltncopyTS,
#else
  .xtrmm_iunucopy = xtrmm_ounucopyTS,
  .xtrmm_iunncopy = xtrmm_ounncopyTS,
  .xtrmm_iutucopy = xtrmm_outucopyTS,
  .xtrmm_iutncopy = xtrmm_outncopyTS,
  .xtrmm_ilnucopy = xtrmm_olnucopyTS,
  .xtrmm_ilnncopy = xtrmm_olnncopyTS,
  .xtrmm_iltucopy = xtrmm_oltucopyTS,
  .xtrmm_iltncopy = xtrmm_oltncopyTS,
#endif
  .xtrmm_ounucopy = xtrmm_ounucopyTS,
  .xtrmm_ounncopy = xtrmm_ounncopyTS,
  .xtrmm_outucopy = xtrmm_outucopyTS,
  .xtrmm_outncopy = xtrmm_outncopyTS,
  .xtrmm_olnucopy = xtrmm_olnucopyTS,
  .xtrmm_olnncopy = xtrmm_olnncopyTS,
  .xtrmm_oltucopy = xtrmm_oltucopyTS,
  .xtrmm_oltncopy = xtrmm_oltncopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_xsymm_dispatch_t openblas_xsymm_dispatchTS = {
#if XGEMM_DEFAULT_UNROLL_M != XGEMM_DEFAULT_UNROLL_N
  .xsymm_iutcopy = xsymm_iutcopyTS,
  .xsymm_iltcopy = xsymm_iltcopyTS,
#else
  .xsymm_iutcopy = xsymm_outcopyTS,
  .xsymm_iltcopy = xsymm_oltcopyTS,
#endif
  .xsymm_outcopy = xsymm_outcopyTS,
  .xsymm_oltcopy = xsymm_oltcopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_xhemm_dispatch_t openblas_xhemm_dispatchTS = {
#if XGEMM_DEFAULT_UNROLL_M != XGEMM_DEFAULT_UNROLL_N
  .xhemm_iutcopy = xhemm_iutcopyTS,
  .xhemm_iltcopy = xhemm_iltcopyTS,
#else
  .xhemm_iutcopy = xhemm_outcopyTS,
  .xhemm_iltcopy = xhemm_oltcopyTS,
#endif
  .xhemm_outcopy = xhemm_outcopyTS,
  .xhemm_oltcopy = xhemm_oltcopyTS,
};
#endif

#ifdef EXPRECISION
const openblas_xgemm3m_dispatch_t openblas_xgemm3m_dispatchTS = {
#if (USE_GEMM3M)
  .xgemm3m_kernel = xgemm3m_kernelTS,
  .xgemm3m_incopyb = xgemm3m_incopybTS,
  .xgemm3m_incopyr = xgemm3m_incopyrTS,
  .xgemm3m_incopyi = xgemm3m_incopyiTS,
  .xgemm3m_itcopyb = xgemm3m_itcopybTS,
  .xgemm3m_itcopyr = xgemm3m_itcopyrTS,
  .xgemm3m_itcopyi = xgemm3m_itcopyiTS,
  .xgemm3m_oncopyb = xgemm3m_oncopybTS,
  .xgemm3m_oncopyr = xgemm3m_oncopyrTS,
  .xgemm3m_oncopyi = xgemm3m_oncopyiTS,
  .xgemm3m_otcopyb = xgemm3m_otcopybTS,
  .xgemm3m_otcopyr = xgemm3m_otcopyrTS,
  .xgemm3m_otcopyi = xgemm3m_otcopyiTS,
#else
  .xgemm3m_kernel = NULL,
  .xgemm3m_incopyb = NULL,
  .xgemm3m_incopyr = NULL,
  .xgemm3m_incopyi = NULL,
  .xgemm3m_itcopyb = NULL,
  .xgemm3m_itcopyr = NULL,
  .xgemm3m_itcopyi = NULL,
  .xgemm3m_oncopyb = NULL,
  .xgemm3m_oncopyr = NULL,
  .xgemm3m_oncopyi = NULL,
  .xgemm3m_otcopyb = NULL,
  .xgemm3m_otcopyr = NULL,
  .xgemm3m_otcopyi = NULL,
#endif
};
#endif

#ifdef EXPRECISION
const openblas_xsymm3m_dispatch_t openblas_xsymm3m_dispatchTS = {
#if (USE_GEMM3M)
  .xsymm3m_iucopyb = xsymm3m_iucopybTS,
  .xsymm3m_ilcopyb = xsymm3m_ilcopybTS,
  .xsymm3m_iucopyr = xsymm3m_iucopyrTS,
  .xsymm3m_ilcopyr = xsymm3m_ilcopyrTS,
  .xsymm3m_iucopyi = xsymm3m_iucopyiTS,
  .xsymm3m_ilcopyi = xsymm3m_ilcopyiTS,
  .xsymm3m_oucopyb = xsymm3m_oucopybTS,
  .xsymm3m_olcopyb = xsymm3m_olcopybTS,
  .xsymm3m_oucopyr = xsymm3m_oucopyrTS,
  .xsymm3m_olcopyr = xsymm3m_olcopyrTS,
  .xsymm3m_oucopyi = xsymm3m_oucopyiTS,
  .xsymm3m_olcopyi = xsymm3m_olcopyiTS,
#else
  .xsymm3m_iucopyb = NULL,
  .xsymm3m_ilcopyb = NULL,
  .xsymm3m_iucopyr = NULL,
  .xsymm3m_ilcopyr = NULL,
  .xsymm3m_iucopyi = NULL,
  .xsymm3m_ilcopyi = NULL,
  .xsymm3m_oucopyb = NULL,
  .xsymm3m_olcopyb = NULL,
  .xsymm3m_oucopyr = NULL,
  .xsymm3m_olcopyr = NULL,
  .xsymm3m_oucopyi = NULL,
  .xsymm3m_olcopyi = NULL,
#endif
};
#endif

#ifdef EXPRECISION
const openblas_xhemm3m_dispatch_t openblas_xhemm3m_dispatchTS = {
#if (USE_GEMM3M)
  .xhemm3m_iucopyb = xhemm3m_iucopybTS,
  .xhemm3m_ilcopyb = xhemm3m_ilcopybTS,
  .xhemm3m_iucopyr = xhemm3m_iucopyrTS,
  .xhemm3m_ilcopyr = xhemm3m_ilcopyrTS,
  .xhemm3m_iucopyi = xhemm3m_iucopyiTS,
  .xhemm3m_ilcopyi = xhemm3m_ilcopyiTS,
  .xhemm3m_oucopyb = xhemm3m_oucopybTS,
  .xhemm3m_olcopyb = xhemm3m_olcopybTS,
  .xhemm3m_oucopyr = xhemm3m_oucopyrTS,
  .xhemm3m_olcopyr = xhemm3m_olcopyrTS,
  .xhemm3m_oucopyi = xhemm3m_oucopyiTS,
  .xhemm3m_olcopyi = xhemm3m_olcopyiTS,
#else
  .xhemm3m_iucopyb = NULL,
  .xhemm3m_ilcopyb = NULL,
  .xhemm3m_iucopyr = NULL,
  .xhemm3m_ilcopyr = NULL,
  .xhemm3m_iucopyi = NULL,
  .xhemm3m_ilcopyi = NULL,
  .xhemm3m_oucopyb = NULL,
  .xhemm3m_olcopyb = NULL,
  .xhemm3m_oucopyr = NULL,
  .xhemm3m_olcopyr = NULL,
  .xhemm3m_oucopyi = NULL,
  .xhemm3m_olcopyi = NULL,
#endif
};
#endif

#ifdef EXPRECISION
const openblas_xneg_dispatch_t openblas_xneg_dispatchTS = {
#ifndef NO_LAPACK
  .xneg_tcopy = xneg_tcopyTS,
#else
  .xneg_tcopy = NULL,
#endif
};
#endif

#ifdef EXPRECISION
const openblas_xlaswp_dispatch_t openblas_xlaswp_dispatchTS = {
#ifndef NO_LAPACK
  .xlaswp_ncopy = xlaswp_ncopyTS,
#else
  .xlaswp_ncopy = NULL,
#endif
};
#endif

#if (BUILD_SINGLE==1)
const openblas_saxpby_dispatch_t openblas_saxpby_dispatchTS = {
  .saxpby_k = saxpby_kTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_daxpby_dispatch_t openblas_daxpby_dispatchTS = {
  .daxpby_k = daxpby_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_caxpby_dispatch_t openblas_caxpby_dispatchTS = {
  .caxpby_k = caxpby_kTS,
};
#endif

#if (BUILD_COMPLEX16==1)
const openblas_zaxpby_dispatch_t openblas_zaxpby_dispatchTS = {
  .zaxpby_k = zaxpby_kTS,
};
#endif

#if (BUILD_SINGLE==1)
const openblas_somatcopy_dispatch_t openblas_somatcopy_dispatchTS = {
  .somatcopy_k_cn = somatcopy_k_cnTS,
  .somatcopy_k_ct = somatcopy_k_ctTS,
  .somatcopy_k_rn = somatcopy_k_rnTS,
  .somatcopy_k_rt = somatcopy_k_rtTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_domatcopy_dispatch_t openblas_domatcopy_dispatchTS = {
  .domatcopy_k_cn = domatcopy_k_cnTS,
  .domatcopy_k_ct = domatcopy_k_ctTS,
  .domatcopy_k_rn = domatcopy_k_rnTS,
  .domatcopy_k_rt = domatcopy_k_rtTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_comatcopy_dispatch_t openblas_comatcopy_dispatchTS = {
  .comatcopy_k_cn = comatcopy_k_cnTS,
  .comatcopy_k_ct = comatcopy_k_ctTS,
  .comatcopy_k_rn = comatcopy_k_rnTS,
  .comatcopy_k_rt = comatcopy_k_rtTS,
  .comatcopy_k_cnc = comatcopy_k_cncTS,
  .comatcopy_k_ctc = comatcopy_k_ctcTS,
  .comatcopy_k_rnc = comatcopy_k_rncTS,
  .comatcopy_k_rtc = comatcopy_k_rtcTS,
};
#endif

#if (BUILD_COMPLEX16==1)
const openblas_zomatcopy_dispatch_t openblas_zomatcopy_dispatchTS = {
  .zomatcopy_k_cn = zomatcopy_k_cnTS,
  .zomatcopy_k_ct = zomatcopy_k_ctTS,
  .zomatcopy_k_rn = zomatcopy_k_rnTS,
  .zomatcopy_k_rt = zomatcopy_k_rtTS,
  .zomatcopy_k_cnc = zomatcopy_k_cncTS,
  .zomatcopy_k_ctc = zomatcopy_k_ctcTS,
  .zomatcopy_k_rnc = zomatcopy_k_rncTS,
  .zomatcopy_k_rtc = zomatcopy_k_rtcTS,
};
#endif

#if (BUILD_SINGLE==1)
const openblas_simatcopy_dispatch_t openblas_simatcopy_dispatchTS = {
  .simatcopy_k_cn = simatcopy_k_cnTS,
  .simatcopy_k_ct = simatcopy_k_ctTS,
  .simatcopy_k_rn = simatcopy_k_rnTS,
  .simatcopy_k_rt = simatcopy_k_rtTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dimatcopy_dispatch_t openblas_dimatcopy_dispatchTS = {
  .dimatcopy_k_cn = dimatcopy_k_cnTS,
  .dimatcopy_k_ct = dimatcopy_k_ctTS,
  .dimatcopy_k_rn = dimatcopy_k_rnTS,
  .dimatcopy_k_rt = dimatcopy_k_rtTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cimatcopy_dispatch_t openblas_cimatcopy_dispatchTS = {
  .cimatcopy_k_cn = cimatcopy_k_cnTS,
  .cimatcopy_k_ct = cimatcopy_k_ctTS,
  .cimatcopy_k_rn = cimatcopy_k_rnTS,
  .cimatcopy_k_rt = cimatcopy_k_rtTS,
  .cimatcopy_k_cnc = cimatcopy_k_cncTS,
  .cimatcopy_k_ctc = cimatcopy_k_ctcTS,
  .cimatcopy_k_rnc = cimatcopy_k_rncTS,
  .cimatcopy_k_rtc = cimatcopy_k_rtcTS,
};
#endif

#if (BUILD_COMPLEX16==1)
const openblas_zimatcopy_dispatch_t openblas_zimatcopy_dispatchTS = {
  .zimatcopy_k_cn = zimatcopy_k_cnTS,
  .zimatcopy_k_ct = zimatcopy_k_ctTS,
  .zimatcopy_k_rn = zimatcopy_k_rnTS,
  .zimatcopy_k_rt = zimatcopy_k_rtTS,
  .zimatcopy_k_cnc = zimatcopy_k_cncTS,
  .zimatcopy_k_ctc = zimatcopy_k_ctcTS,
  .zimatcopy_k_rnc = zimatcopy_k_rncTS,
  .zimatcopy_k_rtc = zimatcopy_k_rtcTS,
};
#endif

#if (BUILD_SINGLE==1)
const openblas_sgeadd_dispatch_t openblas_sgeadd_dispatchTS = {
  .sgeadd_k = sgeadd_kTS,
};
#endif

#if (BUILD_DOUBLE==1)
const openblas_dgeadd_dispatch_t openblas_dgeadd_dispatchTS = {
  .dgeadd_k = dgeadd_kTS,
};
#endif

#if (BUILD_COMPLEX==1)
const openblas_cgeadd_dispatch_t openblas_cgeadd_dispatchTS = {
  .cgeadd_k = cgeadd_kTS,
};
#endif

#if (BUILD_COMPLEX16==1)
const openblas_zgeadd_dispatch_t openblas_zgeadd_dispatchTS = {
  .zgeadd_k = zgeadd_kTS,
};
#endif

#if (ARCH_ARM64)
static void init_parameter(void) {
#if (BUILD_BFLOAT16)
  TABLE_NAME.sbgemm_p = SBGEMM_DEFAULT_P;
  TABLE_NAME.bgemm_p = BGEMM_DEFAULT_P;
#endif
#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE == 1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif

#if (BUILD_BFLOAT16)
  TABLE_NAME.sbgemm_q = SBGEMM_DEFAULT_Q;
  TABLE_NAME.bgemm_q = BGEMM_DEFAULT_Q;
#endif
#if BUILD_SINGLE == 1 || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
#endif
#if BUILD_DOUBLE== 1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
#endif
#if BUILD_COMPLEX== 1
  TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_q = ZGEMM_DEFAULT_Q;
#endif

#if (BUILD_BFLOAT16)
  TABLE_NAME.sbgemm_r = SBGEMM_DEFAULT_R;
  TABLE_NAME.bgemm_r = BGEMM_DEFAULT_R;
#endif
#if BUILD_SINGLE == 1 || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_r = SGEMM_DEFAULT_R;
#endif
#if BUILD_DOUBLE==1  || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_r = DGEMM_DEFAULT_R;
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_r = CGEMM_DEFAULT_R;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_r = ZGEMM_DEFAULT_R;
#endif

#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
  TABLE_NAME.qgemm_q = QGEMM_DEFAULT_Q;
  TABLE_NAME.xgemm_q = XGEMM_DEFAULT_Q;
  TABLE_NAME.qgemm_r = QGEMM_DEFAULT_R;
  TABLE_NAME.xgemm_r = XGEMM_DEFAULT_R;
#endif

#if (USE_GEMM3M)
#ifdef CGEMM3M_DEFAULT_P
  TABLE_NAME.cgemm3m_p = CGEMM3M_DEFAULT_P;
#else
  TABLE_NAME.cgemm3m_p = TABLE_NAME.sgemm_p;
#endif

#ifdef ZGEMM3M_DEFAULT_P
  TABLE_NAME.zgemm3m_p = ZGEMM3M_DEFAULT_P;
#else
  TABLE_NAME.zgemm3m_p = TABLE_NAME.dgemm_p;
#endif

#ifdef CGEMM3M_DEFAULT_Q
  TABLE_NAME.cgemm3m_q = CGEMM3M_DEFAULT_Q;
#else
  TABLE_NAME.cgemm3m_q = TABLE_NAME.sgemm_q;
#endif

#ifdef ZGEMM3M_DEFAULT_Q
  TABLE_NAME.zgemm3m_q = ZGEMM3M_DEFAULT_Q;
#else
  TABLE_NAME.zgemm3m_q = TABLE_NAME.dgemm_q;
#endif

#ifdef CGEMM3M_DEFAULT_R
  TABLE_NAME.cgemm3m_r = CGEMM3M_DEFAULT_R;
#else
  TABLE_NAME.cgemm3m_r = TABLE_NAME.sgemm_r;
#endif

#ifdef ZGEMM3M_DEFAULT_R
  TABLE_NAME.zgemm3m_r = ZGEMM3M_DEFAULT_R;
#else
  TABLE_NAME.zgemm3m_r = TABLE_NAME.dgemm_r;
#endif

#ifdef EXPRECISION
  TABLE_NAME.xgemm3m_p = TABLE_NAME.qgemm_p;
  TABLE_NAME.xgemm3m_q = TABLE_NAME.qgemm_q;
  TABLE_NAME.xgemm3m_r = TABLE_NAME.qgemm_r;
#endif
#endif

}
#else // (ARCH_ARM64)
#if defined(ARCH_MIPS64)
static void init_parameter(void) {
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;

  TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
  TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
  TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;
  TABLE_NAME.zgemm_q = ZGEMM_DEFAULT_Q;

  TABLE_NAME.sgemm_r = SGEMM_DEFAULT_R;
  TABLE_NAME.dgemm_r = 640;
  TABLE_NAME.cgemm_r = CGEMM_DEFAULT_R;
  TABLE_NAME.zgemm_r = ZGEMM_DEFAULT_R;

#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
  TABLE_NAME.qgemm_q = QGEMM_DEFAULT_Q;
  TABLE_NAME.xgemm_q = XGEMM_DEFAULT_Q;
  TABLE_NAME.qgemm_r = QGEMM_DEFAULT_R;
  TABLE_NAME.xgemm_r = XGEMM_DEFAULT_R;
#endif

#if defined(USE_GEMM3M)
#ifdef CGEMM3M_DEFAULT_P
  TABLE_NAME.cgemm3m_p = CGEMM3M_DEFAULT_P;
#else
  TABLE_NAME.cgemm3m_p = TABLE_NAME.sgemm_p;
#endif

#ifdef ZGEMM3M_DEFAULT_P
  TABLE_NAME.zgemm3m_p = ZGEMM3M_DEFAULT_P;
#else
  TABLE_NAME.zgemm3m_p = TABLE_NAME.dgemm_p;
#endif

#ifdef CGEMM3M_DEFAULT_Q
  TABLE_NAME.cgemm3m_q = CGEMM3M_DEFAULT_Q;
#else
  TABLE_NAME.cgemm3m_q = TABLE_NAME.sgemm_q;
#endif

#ifdef ZGEMM3M_DEFAULT_Q
  TABLE_NAME.zgemm3m_q = ZGEMM3M_DEFAULT_Q;
#else
  TABLE_NAME.zgemm3m_q = TABLE_NAME.dgemm_q;
#endif

#ifdef CGEMM3M_DEFAULT_R
  TABLE_NAME.cgemm3m_r = CGEMM3M_DEFAULT_R;
#else
  TABLE_NAME.cgemm3m_r = TABLE_NAME.sgemm_r;
#endif

#ifdef ZGEMM3M_DEFAULT_R
  TABLE_NAME.zgemm3m_r = ZGEMM3M_DEFAULT_R;
#else
  TABLE_NAME.zgemm3m_r = TABLE_NAME.dgemm_r;
#endif

#ifdef EXPRECISION
  TABLE_NAME.xgemm3m_p = TABLE_NAME.qgemm_p;
  TABLE_NAME.xgemm3m_q = TABLE_NAME.qgemm_q;
  TABLE_NAME.xgemm3m_r = TABLE_NAME.qgemm_r;
#endif
#endif
}
#else // (ARCH_MIPS64)
#if (ARCH_LOONGARCH64)
static int get_L3_size() {
  int ret = 0, id = 0x14;
  __asm__ volatile (
    "cpucfg %[ret], %[id]"
    : [ret]"=r"(ret)
    : [id]"r"(id)
    : "memory"
  );
  return ((ret & 0xffff) + 1) * pow(2, ((ret >> 16) & 0xff)) * pow(2, ((ret >> 24) & 0x7f)) / 1024 / 1024; // MB
}
static int get_cpu_prid() {
  int ret = 0, id = 0x0;
  __asm__ volatile (
    "cpucfg %[ret], %[id]"
    : [ret]"=r"(ret)
    : [id]"r"(id)
    : "memory"
  );
  return ret;
}
static void init_parameter(void) {

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_p = SBGEMM_DEFAULT_P;
  TABLE_NAME.bgemm_p = BGEMM_DEFAULT_P;
#endif

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_r = SBGEMM_DEFAULT_R;
  TABLE_NAME.bgemm_r = BGEMM_DEFAULT_R;
#endif

#if defined(LA464)
  int L3_size = get_L3_size();
#ifdef SMP
  if(blas_num_threads == 1){
#endif
    //single thread
    if (L3_size == 32){ // 3C5000 and 3D5000
      TABLE_NAME.sgemm_p = 256;
      TABLE_NAME.sgemm_q = 384;
      TABLE_NAME.sgemm_r = 8192;

      TABLE_NAME.dgemm_p = 112;
      TABLE_NAME.dgemm_q = 289;
      TABLE_NAME.dgemm_r = 4096;

      TABLE_NAME.cgemm_p = 128;
      TABLE_NAME.cgemm_q = 256;
      TABLE_NAME.cgemm_r = 4096;

      TABLE_NAME.zgemm_p = 128;
      TABLE_NAME.zgemm_q = 128;
      TABLE_NAME.zgemm_r = 2048;
    } else { // 3A5000 and 3C5000L
      TABLE_NAME.sgemm_p = 256;
      TABLE_NAME.sgemm_q = 384;
      TABLE_NAME.sgemm_r = 4096;

      TABLE_NAME.dgemm_p = 112;
      TABLE_NAME.dgemm_q = 300;
      TABLE_NAME.dgemm_r = 3024;

      TABLE_NAME.cgemm_p = 128;
      TABLE_NAME.cgemm_q = 256;
      TABLE_NAME.cgemm_r = 2048;

      TABLE_NAME.zgemm_p = 128;
      TABLE_NAME.zgemm_q = 128;
      TABLE_NAME.zgemm_r = 1024;
    }
#ifdef SMP
  }else{
    //multi thread
    if (L3_size == 32){ // 3C5000 and 3D5000
      TABLE_NAME.sgemm_p = 256;
      TABLE_NAME.sgemm_q = 384;
      TABLE_NAME.sgemm_r = 1024;

      TABLE_NAME.dgemm_p = 112;
      TABLE_NAME.dgemm_q = 289;
      TABLE_NAME.dgemm_r = 353;

      TABLE_NAME.cgemm_p = 128;
      TABLE_NAME.cgemm_q = 256;
      TABLE_NAME.cgemm_r = 512;

      TABLE_NAME.zgemm_p = 128;
      TABLE_NAME.zgemm_q = 128;
      TABLE_NAME.zgemm_r = 512;
    } else { // 3A5000 and 3C5000L
      TABLE_NAME.sgemm_p = 256;
      TABLE_NAME.sgemm_q = 384;
      TABLE_NAME.sgemm_r = 2048;

      TABLE_NAME.dgemm_p = 112;
      TABLE_NAME.dgemm_q = 300;
      TABLE_NAME.dgemm_r = 738;

      TABLE_NAME.cgemm_p = 128;
      TABLE_NAME.cgemm_q = 256;
      TABLE_NAME.cgemm_r = 1024;

      TABLE_NAME.zgemm_p = 128;
      TABLE_NAME.zgemm_q = 128;
      TABLE_NAME.zgemm_r = 1024;
    }
  }
#endif
#elif defined(LA264)
  int prid = get_cpu_prid();
  if (prid == 0x0014b020) { //2k3000
        TABLE_NAME.zgemm_p = 128;
        TABLE_NAME.zgemm_q = 176;
        TABLE_NAME.zgemm_r = 360;
  } else {
        TABLE_NAME.zgemm_p = 64;
        TABLE_NAME.zgemm_q = 120;
        TABLE_NAME.zgemm_r = 4096;
  }
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;

  TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
  TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
  TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;

  TABLE_NAME.sgemm_r = SGEMM_DEFAULT_R;
  TABLE_NAME.dgemm_r = DGEMM_DEFAULT_R;
  TABLE_NAME.cgemm_r = CGEMM_DEFAULT_R;
#else
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;

  TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
  TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
  TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;
  TABLE_NAME.zgemm_q = ZGEMM_DEFAULT_Q;

  TABLE_NAME.sgemm_r = SGEMM_DEFAULT_R;
  TABLE_NAME.dgemm_r = DGEMM_DEFAULT_R;
  TABLE_NAME.cgemm_r = CGEMM_DEFAULT_R;
  TABLE_NAME.zgemm_r = ZGEMM_DEFAULT_R;
#endif

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_q = SBGEMM_DEFAULT_Q;
  TABLE_NAME.bgemm_q = BGEMM_DEFAULT_Q;
#endif
}
#else // (ARCH_LOONGARCH64)
#if (ARCH_POWER)
static void init_parameter(void) {

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_p = SBGEMM_DEFAULT_P;
  TABLE_NAME.bgemm_p = BGEMM_DEFAULT_P;
#endif
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_r = SBGEMM_DEFAULT_R;
  TABLE_NAME.bgemm_r = BGEMM_DEFAULT_R;
#endif
  TABLE_NAME.sgemm_r = SGEMM_DEFAULT_R;
  TABLE_NAME.dgemm_r = DGEMM_DEFAULT_R;
  TABLE_NAME.cgemm_r = CGEMM_DEFAULT_R;
  TABLE_NAME.zgemm_r = ZGEMM_DEFAULT_R;


#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_q = SBGEMM_DEFAULT_Q;
  TABLE_NAME.bgemm_q = BGEMM_DEFAULT_Q;
#endif
  TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
  TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
  TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;
  TABLE_NAME.zgemm_q = ZGEMM_DEFAULT_Q;
}
#else //POWER

#if (ARCH_ZARCH)
static void init_parameter(void) {
#ifdef BUILD_BFLOAT16
	TABLE_NAME.sbgemm_p = SBGEMM_DEFAULT_P;
	TABLE_NAME.bgemm_p = BGEMM_DEFAULT_P;
#endif
	TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
	TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
	TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
	TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;

#ifdef BUILD_BFLOAT16
	TABLE_NAME.sbgemm_r = SBGEMM_DEFAULT_R;
	TABLE_NAME.bgemm_r = BGEMM_DEFAULT_R;
#endif
	TABLE_NAME.sgemm_r = SGEMM_DEFAULT_R;
	TABLE_NAME.dgemm_r = DGEMM_DEFAULT_R;
	TABLE_NAME.cgemm_r = CGEMM_DEFAULT_R;
	TABLE_NAME.zgemm_r = ZGEMM_DEFAULT_R;


#ifdef BUILD_BFLOAT16
	TABLE_NAME.sbgemm_q = SBGEMM_DEFAULT_Q;
	TABLE_NAME.bgemm_q = BGEMM_DEFAULT_Q;
#endif
	TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
	TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
	TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;
	TABLE_NAME.zgemm_q = ZGEMM_DEFAULT_Q;
}
#else //ZARCH

#if (ARCH_RISCV64)
static void init_parameter(void) {

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_p = SBGEMM_DEFAULT_P;
  TABLE_NAME.bgemm_p = BGEMM_DEFAULT_P;
#endif
#ifdef BUILD_HFLOAT16
  TABLE_NAME.shgemm_p = SHGEMM_DEFAULT_P;
#endif
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_r = SBGEMM_DEFAULT_R;
  TABLE_NAME.bgemm_r = BGEMM_DEFAULT_R;
#endif
#ifdef BUILD_HFLOAT16
  TABLE_NAME.shgemm_r = SHGEMM_DEFAULT_R;
#endif
  TABLE_NAME.sgemm_r = SGEMM_DEFAULT_R;
  TABLE_NAME.dgemm_r = DGEMM_DEFAULT_R;
  TABLE_NAME.cgemm_r = CGEMM_DEFAULT_R;
  TABLE_NAME.zgemm_r = ZGEMM_DEFAULT_R;


#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_q = SBGEMM_DEFAULT_Q;
  TABLE_NAME.bgemm_q = BGEMM_DEFAULT_Q;
#endif
#ifdef BUILD_HFLOAT16
  TABLE_NAME.shgemm_q = SHGEMM_DEFAULT_Q;
#endif
  TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
  TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
  TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;
  TABLE_NAME.zgemm_q = ZGEMM_DEFAULT_Q;
}
#else //RISCV64

#ifdef ARCH_X86
static int get_l2_size_old(void){
  int i, eax, ebx, ecx, edx, cpuid_level;
  int info[15];

  cpuid(2, &eax, &ebx, &ecx, &edx);

  info[ 0] = BITMASK(eax,  8, 0xff);
  info[ 1] = BITMASK(eax, 16, 0xff);
  info[ 2] = BITMASK(eax, 24, 0xff);

  info[ 3] = BITMASK(ebx,  0, 0xff);
  info[ 4] = BITMASK(ebx,  8, 0xff);
  info[ 5] = BITMASK(ebx, 16, 0xff);
  info[ 6] = BITMASK(ebx, 24, 0xff);

  info[ 7] = BITMASK(ecx,  0, 0xff);
  info[ 8] = BITMASK(ecx,  8, 0xff);
  info[ 9] = BITMASK(ecx, 16, 0xff);
  info[10] = BITMASK(ecx, 24, 0xff);

  info[11] = BITMASK(edx,  0, 0xff);
  info[12] = BITMASK(edx,  8, 0xff);
  info[13] = BITMASK(edx, 16, 0xff);
  info[14] = BITMASK(edx, 24, 0xff);

  for (i = 0; i < 15; i++){

    switch (info[i]){

      /* This table is from http://www.sandpile.org/ia32/cpuid.htm */

    case 0x1a :
      return 96;

    case 0x39 :
    case 0x3b :
    case 0x41 :
    case 0x79 :
    case 0x81 :
      return 128;

    case 0x3a :
      return 192;

    case 0x21 :
    case 0x3c :
    case 0x42 :
    case 0x7a :
    case 0x7e :
    case 0x82 :
      return 256;

    case 0x3d :
      return 384;

    case 0x3e :
    case 0x43 :
    case 0x7b :
    case 0x7f :
    case 0x83 :
    case 0x86 :
      return 512;

    case 0x44 :
    case 0x78 :
    case 0x7c :
    case 0x84 :
    case 0x87 :
      return 1024;

    case 0x45 :
    case 0x7d :
    case 0x85 :
      return 2048;

    case 0x48 :
      return 3184;

    case 0x49 :
      return 4096;

    case 0x4e :
      return 6144;
    }
  }
//  return 0;
fprintf (stderr,"OpenBLAS WARNING - could not determine the L2 cache size on this system, assuming 256k\n");
return 256;
}
#endif

static __inline__ int get_l2_size(void){

  int eax, ebx, ecx, edx, l2;

  l2 = readenv_atoi("OPENBLAS_L2_SIZE");
  if (l2 != 0)
    return l2;

  cpuid(0x80000006, &eax, &ebx, &ecx, &edx);

  l2 = BITMASK(ecx, 16, 0xffff);

#ifndef ARCH_X86
  if (l2 <= 0) {
     fprintf (stderr,"OpenBLAS WARNING - could not determine the L2 cache size on this system, assuming 256k\n");
     return 256;
  }
  return l2;

#else

  if (l2 > 0) return l2;

  return get_l2_size_old();
#endif
}

static __inline__ int get_l3_size(void){

  int eax, ebx, ecx, edx;

  cpuid(0x80000006, &eax, &ebx, &ecx, &edx);

  return BITMASK(edx, 18, 0x3fff) * 512;
}


static void init_parameter(void) {

  int l2 = get_l2_size();

  (void) l2; /* dirty trick to suppress unused variable warning for targets */
             /* where the GEMM unrolling parameters do not depend on l2 */

#ifdef BUILD_BFLOAT16
  TABLE_NAME.sbgemm_p = SBGEMM_DEFAULT_P;
  TABLE_NAME.sbgemm_q = SBGEMM_DEFAULT_Q;
  TABLE_NAME.bgemm_p = BGEMM_DEFAULT_P;
  TABLE_NAME.bgemm_q = BGEMM_DEFAULT_Q;
#endif
#ifdef BUILD_HFLOAT16
  TABLE_NAME.shgemm_p = SHGEMM_DEFAULT_P;
  TABLE_NAME.shgemm_q = SHGEMM_DEFAULT_Q;
#endif
#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_q = SGEMM_DEFAULT_Q;
#endif
#if  (BUILD_DOUBLE==1) || (BUILD_COMPLEX16)
  TABLE_NAME.dgemm_q = DGEMM_DEFAULT_Q;
#endif
#if BUILD_COMPLEX == 1
  TABLE_NAME.cgemm_q = CGEMM_DEFAULT_Q;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_q = ZGEMM_DEFAULT_Q;
#endif

#if BUILD_COMPLEX == 1
#ifdef CGEMM3M_DEFAULT_Q
  TABLE_NAME.cgemm3m_q = CGEMM3M_DEFAULT_Q;
#else
  TABLE_NAME.cgemm3m_q = SGEMM_DEFAULT_Q;
#endif
#endif

#if BUILD_COMPLEX16 == 1
#ifdef ZGEMM3M_DEFAULT_Q
  TABLE_NAME.zgemm3m_q = ZGEMM3M_DEFAULT_Q;
#else
  TABLE_NAME.zgemm3m_q = DGEMM_DEFAULT_Q;
#endif
#endif

#ifdef EXPRECISION
  TABLE_NAME.qgemm_q = QGEMM_DEFAULT_Q;
  TABLE_NAME.xgemm_q = XGEMM_DEFAULT_Q;
  TABLE_NAME.xgemm3m_q = QGEMM_DEFAULT_Q;
#endif

#if defined(CORE_KATMAI)  || defined(CORE_COPPERMINE) || defined(CORE_BANIAS) || defined(CORE_YONAH) || defined(CORE_ATHLON)

#ifdef DEBUG
  fprintf(stderr, "Katmai, Coppermine, Banias, Athlon\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p =  64 * (l2 >> 7);
#endif
#if BUILD_DOUBLE == 1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p =  32 * (l2 >> 7);
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p =  32 * (l2 >> 7);
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p =  16 * (l2 >> 7);
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  16 * (l2 >> 7);
  TABLE_NAME.xgemm_p =   8 * (l2 >> 7);
#endif
#endif

#ifdef CORE_NORTHWOOD

#ifdef DEBUG
  fprintf(stderr, "Northwood\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p =  96 * (l2 >> 7);
#endif
#if BUILD_DOUBLE == 1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p =  48 * (l2 >> 7);
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p =  48 * (l2 >> 7);
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p =  24 * (l2 >> 7);
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  24 * (l2 >> 7);
  TABLE_NAME.xgemm_p =  12 * (l2 >> 7);
#endif
#endif

#ifdef ATOM

#ifdef DEBUG
  fprintf(stderr, "Atom\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = 256;
#endif
#if BUILD_DOUBLE ==1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = 128;
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p = 128;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p =  64;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  64;
  TABLE_NAME.xgemm_p =  32;
#endif
#endif

#ifdef CORE_PRESCOTT

#ifdef DEBUG
  fprintf(stderr, "Prescott\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p =  56 * (l2 >> 7);
#endif
#if BUILD_DOUBLE ==1  || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p =  28 * (l2 >> 7);
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p =  28 * (l2 >> 7);
#endif
#if BUILD_COMPLEX16 == 1
  TABLE_NAME.zgemm_p =  14 * (l2 >> 7);
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  14 * (l2 >> 7);
  TABLE_NAME.xgemm_p =   7 * (l2 >> 7);
#endif
#endif

#ifdef CORE2

#ifdef DEBUG
  fprintf(stderr, "Core2\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p =  92 * (l2 >> 9) + 8;
#endif
#if BUILD_DOUBLE==1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p =  46 * (l2 >> 9) + 8;
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p =  46 * (l2 >> 9) + 4;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p =  23 * (l2 >> 9) + 4;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  92 * (l2 >> 9) + 8;
  TABLE_NAME.xgemm_p =  46 * (l2 >> 9) + 4;
#endif
#endif

#ifdef PENRYN

#ifdef DEBUG
  fprintf(stderr, "Penryn\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p =  42 * (l2 >> 9) + 8;
#endif
#if BUILD_DOUBLE == 1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p =  42 * (l2 >> 9) + 8;
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p =  21 * (l2 >> 9) + 4;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p =  21 * (l2 >> 9) + 4;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  42 * (l2 >> 9) + 8;
  TABLE_NAME.xgemm_p =  21 * (l2 >> 9) + 4;
#endif
#endif

#ifdef DUNNINGTON

#ifdef DEBUG
  fprintf(stderr, "Dunnington\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p =  42 * (l2 >> 9) + 8;
#endif
#if BUILD_DOUBLE ==1 || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p =  42 * (l2 >> 9) + 8;
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p =  21 * (l2 >> 9) + 4;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p =  21 * (l2 >> 9) + 4;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  42 * (l2 >> 9) + 8;
  TABLE_NAME.xgemm_p =  21 * (l2 >> 9) + 4;
#endif
#endif


#ifdef NEHALEM

#ifdef DEBUG
  fprintf(stderr, "Nehalem\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#ifdef SANDYBRIDGE

#ifdef DEBUG
  fprintf(stderr, "Sandybridge\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#ifdef HASWELL

#ifdef DEBUG
  fprintf(stderr, "Haswell\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)

#ifdef DEBUG
  fprintf(stderr, "SkylakeX\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif


#ifdef OPTERON

#ifdef DEBUG
  fprintf(stderr, "Opteron\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = 224 +  56 * (l2 >> 7);
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = 112 +  28 * (l2 >> 7);
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = 112 +  28 * (l2 >> 7);
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p =  56 +  14 * (l2 >> 7);
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p =  56 +  14 * (l2 >> 7);
  TABLE_NAME.xgemm_p =  28 +   7 * (l2 >> 7);
#endif
#endif

#ifdef BARCELONA

#ifdef DEBUG
  fprintf(stderr, "Barcelona\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#ifdef BOBCAT

#ifdef DEBUG
  fprintf(stderr, "Bobcate\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#ifdef BULLDOZER

#ifdef DEBUG
  fprintf(stderr, "Bulldozer\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#ifdef EXCAVATOR

#ifdef DEBUG
  fprintf(stderr, "Excavator\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif


#ifdef PILEDRIVER

#ifdef DEBUG
  fprintf(stderr, "Piledriver\n");
#endif

#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#ifdef STEAMROLLER

#ifdef DEBUG
  fprintf(stderr, "Steamroller\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if BUILD_DOUBLE || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif

#ifdef ZEN

#ifdef DEBUG
  fprintf(stderr, "Zen\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if BUILD_COMPLEX16
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif
#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif
#endif


#ifdef NANO

#ifdef DEBUG
  fprintf(stderr, "NANO\n");
#endif

#if  (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
  TABLE_NAME.sgemm_p = SGEMM_DEFAULT_P;
#endif
#if  (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
  TABLE_NAME.dgemm_p = DGEMM_DEFAULT_P;
#endif
#if (BUILD_COMPLEX==1)
  TABLE_NAME.cgemm_p = CGEMM_DEFAULT_P;
#endif
#if (BUILD_COMPLEX16==1)
  TABLE_NAME.zgemm_p = ZGEMM_DEFAULT_P;
#endif


#ifdef EXPRECISION
  TABLE_NAME.qgemm_p = QGEMM_DEFAULT_P;
  TABLE_NAME.xgemm_p = XGEMM_DEFAULT_P;
#endif

#endif

#ifdef SAPPHIRERAPIDS
#if (BUILD_BFLOAT16 == 1)
  TABLE_NAME.need_amxtile_permission = 1;
#endif
#endif

#if BUILD_COMPLEX==1
#ifdef CGEMM3M_DEFAULT_P
  TABLE_NAME.cgemm3m_p = CGEMM3M_DEFAULT_P;
#else
  TABLE_NAME.cgemm3m_p = TABLE_NAME.sgemm_p;
#endif
#endif

#if BUILD_COMPLEX16==1
#ifdef ZGEMM3M_DEFAULT_P
  TABLE_NAME.zgemm3m_p = ZGEMM3M_DEFAULT_P;
#else
  TABLE_NAME.zgemm3m_p = TABLE_NAME.dgemm_p;
#endif
#endif

#ifdef EXPRECISION
  TABLE_NAME.xgemm3m_p = TABLE_NAME.qgemm_p;
#endif
	
#ifndef NO_AVX512
{
    int l3_kb = get_l3_size();
    int l2_kb = get_l2_size();
    unsigned int eax, ebx, ecx, edx;
    unsigned int cpuid7_eax, cpuid7_ebx, cpuid7_ecx, cpuid7_edx;

    cpuid(0, &eax, &ebx, &ecx, &edx);

    if ((ebx == 0x68747541) && (l3_kb > 0) && (l3_kb % 32768 == 0) && (l2_kb == 1024)) { //Auth AMD
      if (strcmp(gotoblas_corename(), "cooperlake") == 0 || strcmp(gotoblas_corename(), "skylakex") == 0 || strcmp(gotoblas_corename(), "sapphirerapids") == 0) {

        cpuid(7, &cpuid7_eax, &cpuid7_ebx, &cpuid7_ecx, &cpuid7_edx);
        
        if (cpuid7_ebx & (1 << 16)) { // avx512 - Zen 4, 5
#if BUILD_SINGLE == 1
            TABLE_NAME.sgemm_p = 384;
            TABLE_NAME.sgemm_q = 512;
#endif
#if BUILD_DOUBLE == 1
            TABLE_NAME.dgemm_p = 512;
            TABLE_NAME.dgemm_q = 512;
#endif
#if BUILD_COMPLEX == 1
            TABLE_NAME.cgemm_p = 160;
            TABLE_NAME.cgemm_q = 480;
#endif
#if BUILD_COMPLEX16 == 1
            TABLE_NAME.zgemm_p = 176;
            TABLE_NAME.zgemm_q = 256;
#endif
        }
    }
  }
}
#endif

#if BUILD_SINGLE == 1
  TABLE_NAME.sgemm_p = ((TABLE_NAME.sgemm_p + SGEMM_DEFAULT_UNROLL_M - 1)/SGEMM_DEFAULT_UNROLL_M) * SGEMM_DEFAULT_UNROLL_M;
#endif
#if BUILD_DOUBLE== 1
  TABLE_NAME.dgemm_p = ((TABLE_NAME.dgemm_p + DGEMM_DEFAULT_UNROLL_M - 1)/DGEMM_DEFAULT_UNROLL_M) * DGEMM_DEFAULT_UNROLL_M;
#endif
#if BUILD_COMPLEX==1
  TABLE_NAME.cgemm_p = ((TABLE_NAME.cgemm_p + CGEMM_DEFAULT_UNROLL_M - 1)/CGEMM_DEFAULT_UNROLL_M) * CGEMM_DEFAULT_UNROLL_M;
#endif
#if BUILD_COMPLEX16==1
  TABLE_NAME.zgemm_p = ((TABLE_NAME.zgemm_p + ZGEMM_DEFAULT_UNROLL_M - 1)/ZGEMM_DEFAULT_UNROLL_M) * ZGEMM_DEFAULT_UNROLL_M;
#endif

#if BUILD_COMPLEX==1
#ifdef CGEMM3M_DEFAULT_UNROLL_M
  TABLE_NAME.cgemm3m_p = ((TABLE_NAME.cgemm3m_p + CGEMM3M_DEFAULT_UNROLL_M - 1)/CGEMM3M_DEFAULT_UNROLL_M) * CGEMM3M_DEFAULT_UNROLL_M;
#else
  TABLE_NAME.cgemm3m_p = ((TABLE_NAME.cgemm3m_p + SGEMM_DEFAULT_UNROLL_M - 1)/SGEMM_DEFAULT_UNROLL_M) * SGEMM_DEFAULT_UNROLL_M;
#endif
#endif

#if BUILD_COMPLEX16==1
#ifdef ZGEMM3M_DEFAULT_UNROLL_M
  TABLE_NAME.zgemm3m_p = ((TABLE_NAME.zgemm3m_p + ZGEMM3M_DEFAULT_UNROLL_M - 1)/ZGEMM3M_DEFAULT_UNROLL_M) * ZGEMM3M_DEFAULT_UNROLL_M;
#else
  TABLE_NAME.zgemm3m_p = ((TABLE_NAME.zgemm3m_p + DGEMM_DEFAULT_UNROLL_M - 1)/DGEMM_DEFAULT_UNROLL_M) * DGEMM_DEFAULT_UNROLL_M;
#endif
#endif

#ifdef QUAD_PRECISION
  TABLE_NAME.qgemm_p = ((TABLE_NAME.qgemm_p + QGEMM_DEFAULT_UNROLL_M - 1)/QGEMM_DEFAULT_UNROLL_M) * QGEMM_DEFAULT_UNROLL_M;
  TABLE_NAME.xgemm_p = ((TABLE_NAME.xgemm_p + XGEMM_DEFAULT_UNROLL_M - 1)/XGEMM_DEFAULT_UNROLL_M) * XGEMM_DEFAULT_UNROLL_M;
  TABLE_NAME.xgemm3m_p = ((TABLE_NAME.xgemm3m_p + QGEMM_DEFAULT_UNROLL_M - 1)/QGEMM_DEFAULT_UNROLL_M) * QGEMM_DEFAULT_UNROLL_M;
#endif

#ifdef DEBUG
  fprintf(stderr, "L2 = %8d DGEMM_P  .. %d\n", l2, TABLE_NAME.dgemm_p);
#endif

#if BUILD_BFLOAT16==1
  TABLE_NAME.sbgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.sbgemm_p * TABLE_NAME.sbgemm_q *  4 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.sbgemm_q *  4) - 15) & ~15);
  TABLE_NAME.bgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.bgemm_p * TABLE_NAME.bgemm_q *  4 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.bgemm_q *  4) - 15) & ~15);
#endif

#if BUILD_HFLOAT16==1
  TABLE_NAME.shgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.shgemm_p * TABLE_NAME.shgemm_q *  4 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.shgemm_q *  4) - 15) & ~15);
#endif

#if BUILD_SINGLE==1
  TABLE_NAME.sgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.sgemm_p * TABLE_NAME.sgemm_q *  4 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.sgemm_q *  4) - 15) & ~15);
#endif

#if BUILD_DOUBLE==1
  TABLE_NAME.dgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.dgemm_p * TABLE_NAME.dgemm_q *  8 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.dgemm_q *  8) - 15) & ~15);
#endif

#ifdef EXPRECISION
  TABLE_NAME.qgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.qgemm_p * TABLE_NAME.qgemm_q * 16 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.qgemm_q * 16) - 15) & ~15);
#endif

#if BUILD_COMPLEX ==1
  TABLE_NAME.cgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.cgemm_p * TABLE_NAME.cgemm_q *  8 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.cgemm_q *  8) - 15) & ~15);
#endif

#if BUILD_COMPLEX16 ==1
  TABLE_NAME.zgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.zgemm_p * TABLE_NAME.zgemm_q * 16 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.zgemm_q * 16) - 15) & ~15);
#endif

#if BUILD_COMPLEX == 1
  TABLE_NAME.cgemm3m_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.cgemm3m_p * TABLE_NAME.cgemm3m_q *  8 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.cgemm3m_q *  8) - 15) & ~15);
#endif

#if BUILD_COMPLEX16 == 1
  TABLE_NAME.zgemm3m_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.zgemm3m_p * TABLE_NAME.zgemm3m_q * 16 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
			       ) / (TABLE_NAME.zgemm3m_q * 16) - 15) & ~15);
#endif



#ifdef EXPRECISION
  TABLE_NAME.xgemm_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.xgemm_p * TABLE_NAME.xgemm_q * 32 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
		       ) / (TABLE_NAME.xgemm_q * 32) - 15) & ~15);

  TABLE_NAME.xgemm3m_r = (((BUFFER_SIZE -
			       ((TABLE_NAME.xgemm3m_p * TABLE_NAME.xgemm3m_q * 32 + TABLE_NAME.offsetA
				 + TABLE_NAME.align) & ~TABLE_NAME.align)
		       ) / (TABLE_NAME.xgemm3m_q * 32) - 15) & ~15);

#endif



}
#endif //RISCV64
#endif //POWER
#endif //ZARCH
#endif //(ARCH_LOONGARCH64)
#endif //(ARCH_MIPS64)
#endif //(ARCH_ARM64)
