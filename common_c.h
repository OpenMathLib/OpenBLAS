/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright 2025-2026 The OpenBLAS Project.                         */
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

#ifndef COMMON_C_H
#define COMMON_C_H

#ifndef DYNAMIC_ARCH

#define	CAMAX_K			camax_k
#define	CAMIN_K			camin_k
#define	CMAX_K			cmax_k
#define	CMIN_K			cmin_k
#define	ICAMAX_K		icamax_k
#define	ICAMIN_K		icamin_k
#define	ICMAX_K			icmax_k
#define	ICMIN_K			icmin_k
#define	CASUM_K			casum_k
#define	CAXPYU_K		caxpy_k
#define	CAXPYC_K		caxpyc_k
#define	CCOPY_K			ccopy_k
#define	CDOTU_K			cdotu_k
#define	CDOTC_K			cdotc_k
#define	CNRM2_K			cnrm2_k
#define	CSCAL_K			cscal_k
#define	CSUM_K			csum_k
#define	CSWAP_K			cswap_k
#define	CROT_K			csrot_k

#define	CGEMV_N			cgemv_n
#define	CGEMV_T			cgemv_t
#define	CGEMV_R			cgemv_r
#define	CGEMV_C			cgemv_c
#define	CGEMV_O			cgemv_o
#define	CGEMV_U			cgemv_u
#define	CGEMV_S			cgemv_s
#define	CGEMV_D			cgemv_d

#define	CGERU_K			cgeru_k
#define	CGERC_K			cgerc_k
#define	CGERV_K			cgerv_k
#define	CGERD_K			cgerd_k

#define CSYMV_U			csymv_U
#define CSYMV_L			csymv_L
#define CHEMV_U			chemv_U
#define CHEMV_L			chemv_L
#define CHEMV_V			chemv_V
#define CHEMV_M			chemv_M

#define CSYMV_THREAD_U		csymv_thread_U
#define CSYMV_THREAD_L		csymv_thread_L
#define CHEMV_THREAD_U		chemv_thread_U
#define CHEMV_THREAD_L		chemv_thread_L
#define CHEMV_THREAD_V		chemv_thread_V
#define CHEMV_THREAD_M		chemv_thread_M

#define	CGEMM_ONCOPY		cgemm_oncopy
#define	CGEMM_OTCOPY		cgemm_otcopy

#if CGEMM_DEFAULT_UNROLL_M == CGEMM_DEFAULT_UNROLL_N
#define	CGEMM_INCOPY		cgemm_oncopy
#define	CGEMM_ITCOPY		cgemm_otcopy
#else
#define	CGEMM_INCOPY		cgemm_incopy
#define	CGEMM_ITCOPY		cgemm_itcopy
#endif

#define	CSYMM_INCOPY		    csymm_incopy
#define	CSYMM_ITCOPY	    	csymm_itcopy
#define	CTRMM_INCOPY		ctrmm_incopy
#define	CTRMM_ITCOPY		ctrmm_itcopy

#define	CTRMM_OUNUCOPY		ctrmm_ounucopy
#define	CTRMM_OUNNCOPY		ctrmm_ounncopy
#define	CTRMM_OUTUCOPY		ctrmm_outucopy
#define	CTRMM_OUTNCOPY		ctrmm_outncopy
#define	CTRMM_OLNUCOPY		ctrmm_olnucopy
#define	CTRMM_OLNNCOPY		ctrmm_olnncopy
#define	CTRMM_OLTUCOPY		ctrmm_oltucopy
#define	CTRMM_OLTNCOPY		ctrmm_oltncopy

#define	CTRSM_OUNUCOPY		ctrsm_ounucopy
#define	CTRSM_OUNNCOPY		ctrsm_ounncopy
#define	CTRSM_OUTUCOPY		ctrsm_outucopy
#define	CTRSM_OUTNCOPY		ctrsm_outncopy
#define	CTRSM_OLNUCOPY		ctrsm_olnucopy
#define	CTRSM_OLNNCOPY		ctrsm_olnncopy
#define	CTRSM_OLTUCOPY		ctrsm_oltucopy
#define	CTRSM_OLTNCOPY		ctrsm_oltncopy

#if CGEMM_DEFAULT_UNROLL_M == CGEMM_DEFAULT_UNROLL_N
#define	CTRMM_IUNUCOPY		ctrmm_ounucopy
#define	CTRMM_IUNNCOPY		ctrmm_ounncopy
#define	CTRMM_IUTUCOPY		ctrmm_outucopy
#define	CTRMM_IUTNCOPY		ctrmm_outncopy
#define	CTRMM_ILNUCOPY		ctrmm_olnucopy
#define	CTRMM_ILNNCOPY		ctrmm_olnncopy
#define	CTRMM_ILTUCOPY		ctrmm_oltucopy
#define	CTRMM_ILTNCOPY		ctrmm_oltncopy

#define	CTRSM_IUNUCOPY		ctrsm_ounucopy
#define	CTRSM_IUNNCOPY		ctrsm_ounncopy
#define	CTRSM_IUTUCOPY		ctrsm_outucopy
#define	CTRSM_IUTNCOPY		ctrsm_outncopy
#define	CTRSM_ILNUCOPY		ctrsm_olnucopy
#define	CTRSM_ILNNCOPY		ctrsm_olnncopy
#define	CTRSM_ILTUCOPY		ctrsm_oltucopy
#define	CTRSM_ILTNCOPY		ctrsm_oltncopy
#else
#define	CTRMM_IUNUCOPY		ctrmm_iunucopy
#define	CTRMM_IUNNCOPY		ctrmm_iunncopy
#define	CTRMM_IUTUCOPY		ctrmm_iutucopy
#define	CTRMM_IUTNCOPY		ctrmm_iutncopy
#define	CTRMM_ILNUCOPY		ctrmm_ilnucopy
#define	CTRMM_ILNNCOPY		ctrmm_ilnncopy
#define	CTRMM_ILTUCOPY		ctrmm_iltucopy
#define	CTRMM_ILTNCOPY		ctrmm_iltncopy

#define	CTRSM_IUNUCOPY		ctrsm_iunucopy
#define	CTRSM_IUNNCOPY		ctrsm_iunncopy
#define	CTRSM_IUTUCOPY		ctrsm_iutucopy
#define	CTRSM_IUTNCOPY		ctrsm_iutncopy
#define	CTRSM_ILNUCOPY		ctrsm_ilnucopy
#define	CTRSM_ILNNCOPY		ctrsm_ilnncopy
#define	CTRSM_ILTUCOPY		ctrsm_iltucopy
#define	CTRSM_ILTNCOPY		ctrsm_iltncopy
#endif

#define	CGEMM_BETA		cgemm_beta
#define	SME_CGEMM_KERNEL	sme_cgemm_kernel

#define	CGEMM_KERNEL_N		cgemm_kernel_n
#define	CGEMM_KERNEL_L		cgemm_kernel_l
#define	CGEMM_KERNEL_R		cgemm_kernel_r
#define	CGEMM_KERNEL_B		cgemm_kernel_b

#define	CSYMM_KERNEL_N		csymm_kernel_n
#define	CSYMM_KERNEL_L		csymm_kernel_l
#define	CSYMM_KERNEL_R		csymm_kernel_r
#define	CSYMM_KERNEL_B		csymm_kernel_b
#define	CTRMM_GEMM_KERNEL_N	ctrmm_gemm_kernel_n
#define	CTRMM_GEMM_KERNEL_L	ctrmm_gemm_kernel_l
#define	CTRMM_GEMM_KERNEL_R	ctrmm_gemm_kernel_r
#define	CTRMM_GEMM_KERNEL_B	ctrmm_gemm_kernel_b

#define	CTRMM_KERNEL_LN		ctrmm_kernel_LN
#define	CTRMM_KERNEL_LT		ctrmm_kernel_LT
#define	CTRMM_KERNEL_LR		ctrmm_kernel_LR
#define	CTRMM_KERNEL_LC		ctrmm_kernel_LC
#define	CTRMM_KERNEL_RN		ctrmm_kernel_RN
#define	CTRMM_KERNEL_RT		ctrmm_kernel_RT
#define	CTRMM_KERNEL_RR		ctrmm_kernel_RR
#define	CTRMM_KERNEL_RC		ctrmm_kernel_RC

#define	CTRSM_KERNEL_LN		ctrsm_kernel_LN
#define	CTRSM_KERNEL_LT		ctrsm_kernel_LT
#define	CTRSM_KERNEL_LR		ctrsm_kernel_LR
#define	CTRSM_KERNEL_LC		ctrsm_kernel_LC
#define	CTRSM_KERNEL_RN		ctrsm_kernel_RN
#define	CTRSM_KERNEL_RT		ctrsm_kernel_RT
#define	CTRSM_KERNEL_RR		ctrsm_kernel_RR
#define	CTRSM_KERNEL_RC		ctrsm_kernel_RC

#define	CSYMM_OUTCOPY		csymm_outcopy
#define	CSYMM_OLTCOPY		csymm_oltcopy
#if CGEMM_DEFAULT_UNROLL_M == CGEMM_DEFAULT_UNROLL_N
#define	CSYMM_IUTCOPY		csymm_outcopy
#define	CSYMM_ILTCOPY		csymm_oltcopy
#else
#define	CSYMM_IUTCOPY		csymm_iutcopy
#define	CSYMM_ILTCOPY		csymm_iltcopy
#endif

#define	CHEMM_OUTCOPY		chemm_outcopy
#define	CHEMM_OLTCOPY		chemm_oltcopy
#if CGEMM_DEFAULT_UNROLL_M == CGEMM_DEFAULT_UNROLL_N
#define	CHEMM_IUTCOPY		chemm_outcopy
#define	CHEMM_ILTCOPY		chemm_oltcopy
#else
#define	CHEMM_IUTCOPY		chemm_iutcopy
#define	CHEMM_ILTCOPY		chemm_iltcopy
#endif

#define	CGEMM3M_ONCOPYB		cgemm3m_oncopyb
#define	CGEMM3M_ONCOPYR		cgemm3m_oncopyr
#define	CGEMM3M_ONCOPYI		cgemm3m_oncopyi
#define	CGEMM3M_OTCOPYB		cgemm3m_otcopyb
#define	CGEMM3M_OTCOPYR		cgemm3m_otcopyr
#define	CGEMM3M_OTCOPYI		cgemm3m_otcopyi

#define	CGEMM3M_INCOPYB		cgemm3m_incopyb
#define	CGEMM3M_INCOPYR		cgemm3m_incopyr
#define	CGEMM3M_INCOPYI		cgemm3m_incopyi
#define	CGEMM3M_ITCOPYB		cgemm3m_itcopyb
#define	CGEMM3M_ITCOPYR		cgemm3m_itcopyr
#define	CGEMM3M_ITCOPYI		cgemm3m_itcopyi

#define	CSYMM3M_ILCOPYB		csymm3m_ilcopyb
#define	CSYMM3M_IUCOPYB		csymm3m_iucopyb
#define	CSYMM3M_ILCOPYR		csymm3m_ilcopyr
#define	CSYMM3M_IUCOPYR		csymm3m_iucopyr
#define	CSYMM3M_ILCOPYI		csymm3m_ilcopyi
#define	CSYMM3M_IUCOPYI		csymm3m_iucopyi

#define	CSYMM3M_OLCOPYB		csymm3m_olcopyb
#define	CSYMM3M_OUCOPYB		csymm3m_oucopyb
#define	CSYMM3M_OLCOPYR		csymm3m_olcopyr
#define	CSYMM3M_OUCOPYR		csymm3m_oucopyr
#define	CSYMM3M_OLCOPYI		csymm3m_olcopyi
#define	CSYMM3M_OUCOPYI		csymm3m_oucopyi

#define	CHEMM3M_ILCOPYB		chemm3m_ilcopyb
#define	CHEMM3M_IUCOPYB		chemm3m_iucopyb
#define	CHEMM3M_ILCOPYR		chemm3m_ilcopyr
#define	CHEMM3M_IUCOPYR		chemm3m_iucopyr
#define	CHEMM3M_ILCOPYI		chemm3m_ilcopyi
#define	CHEMM3M_IUCOPYI		chemm3m_iucopyi

#define	CHEMM3M_OLCOPYB		chemm3m_olcopyb
#define	CHEMM3M_OUCOPYB		chemm3m_oucopyb
#define	CHEMM3M_OLCOPYR		chemm3m_olcopyr
#define	CHEMM3M_OUCOPYR		chemm3m_oucopyr
#define	CHEMM3M_OLCOPYI		chemm3m_olcopyi
#define	CHEMM3M_OUCOPYI		chemm3m_oucopyi

#define	CGEMM3M_KERNEL		cgemm3m_kernel

#define CNEG_TCOPY		cneg_tcopy
#define CLASWP_NCOPY		claswp_ncopy

#define CAXPBY_K                caxpby_k

#define COMATCOPY_K_CN          comatcopy_k_cn
#define COMATCOPY_K_RN          comatcopy_k_rn
#define COMATCOPY_K_CT          comatcopy_k_ct
#define COMATCOPY_K_RT          comatcopy_k_rt
#define COMATCOPY_K_CNC         comatcopy_k_cnc
#define COMATCOPY_K_RNC         comatcopy_k_rnc
#define COMATCOPY_K_CTC         comatcopy_k_ctc
#define COMATCOPY_K_RTC         comatcopy_k_rtc

#define CIMATCOPY_K_CN          cimatcopy_k_cn
#define CIMATCOPY_K_RN          cimatcopy_k_rn
#define CIMATCOPY_K_CT          cimatcopy_k_ct
#define CIMATCOPY_K_RT          cimatcopy_k_rt
#define CIMATCOPY_K_CNC         cimatcopy_k_cnc
#define CIMATCOPY_K_RNC         cimatcopy_k_rnc
#define CIMATCOPY_K_CTC         cimatcopy_k_ctc
#define CIMATCOPY_K_RTC         cimatcopy_k_rtc

#define CGEADD_K                cgeadd_k 

#define CGEMM_SMALL_MATRIX_PERMIT	cgemm_small_matrix_permit

#else

#define	CAMAX_K			OPENBLAS_DISPATCH(camax) -> camax_k
#define	CAMIN_K			OPENBLAS_DISPATCH(camin) -> camin_k
#define	CMAX_K			gotoblas -> cmax_k
#define	CMIN_K			gotoblas -> cmin_k
#define	ICAMAX_K		OPENBLAS_DISPATCH(icamax) -> icamax_k
#define	ICAMIN_K		OPENBLAS_DISPATCH(icamin) -> icamin_k
#define	ICMAX_K			gotoblas -> icmax_k
#define	ICMIN_K			gotoblas -> icmin_k
#define	CASUM_K			OPENBLAS_DISPATCH(casum) -> casum_k
#define	CAXPYU_K		OPENBLAS_DISPATCH(caxpy) -> caxpy_k
#define	CAXPYC_K		OPENBLAS_DISPATCH(caxpyc) -> caxpyc_k
#define	CCOPY_K			OPENBLAS_DISPATCH(ccopy) -> ccopy_k
#define	CDOTU_K			OPENBLAS_DISPATCH(cdotu) -> cdotu_k
#define	CDOTC_K			OPENBLAS_DISPATCH(cdotc) -> cdotc_k
#define	CNRM2_K			OPENBLAS_DISPATCH(cnrm2) -> cnrm2_k
#define	CSCAL_K			OPENBLAS_DISPATCH(cscal) -> cscal_k
#define	CSUM_K			OPENBLAS_DISPATCH(csum) -> csum_k
#define	CSWAP_K			OPENBLAS_DISPATCH(cswap) -> cswap_k
#define	CROT_K			OPENBLAS_DISPATCH(csrot) -> csrot_k

#define	CGEMV_N			OPENBLAS_DISPATCH(cgemv) -> cgemv_n
#define	CGEMV_T			OPENBLAS_DISPATCH(cgemv) -> cgemv_t
#define	CGEMV_R			OPENBLAS_DISPATCH(cgemv) -> cgemv_r
#define	CGEMV_C			OPENBLAS_DISPATCH(cgemv) -> cgemv_c
#define	CGEMV_O			OPENBLAS_DISPATCH(cgemv) -> cgemv_o
#define	CGEMV_U			OPENBLAS_DISPATCH(cgemv) -> cgemv_u
#define	CGEMV_S			OPENBLAS_DISPATCH(cgemv) -> cgemv_s
#define	CGEMV_D			OPENBLAS_DISPATCH(cgemv) -> cgemv_d

#define	CGERU_K			OPENBLAS_DISPATCH(cgeru) -> cgeru_k
#define	CGERC_K			OPENBLAS_DISPATCH(cgerc) -> cgerc_k
#define	CGERV_K			OPENBLAS_DISPATCH(cgerv) -> cgerv_k
#define	CGERD_K			OPENBLAS_DISPATCH(cgerd) -> cgerd_k

#define CSYMV_U			OPENBLAS_DISPATCH(csymv) -> csymv_U
#define CSYMV_L			OPENBLAS_DISPATCH(csymv) -> csymv_L
#define CHEMV_U			OPENBLAS_DISPATCH(chemv) -> chemv_U
#define CHEMV_L			OPENBLAS_DISPATCH(chemv) -> chemv_L
#define CHEMV_V			OPENBLAS_DISPATCH(chemv) -> chemv_V
#define CHEMV_M			OPENBLAS_DISPATCH(chemv) -> chemv_M

#define CSYMV_THREAD_U		csymv_thread_U
#define CSYMV_THREAD_L		csymv_thread_L
#define CHEMV_THREAD_U		chemv_thread_U
#define CHEMV_THREAD_L		chemv_thread_L
#define CHEMV_THREAD_V		chemv_thread_V
#define CHEMV_THREAD_M		chemv_thread_M

#define	CGEMM_ONCOPY		OPENBLAS_DISPATCH(cgemm) -> cgemm_oncopy
#define	CGEMM_OTCOPY		OPENBLAS_DISPATCH(cgemm) -> cgemm_otcopy
#define	CGEMM_INCOPY		OPENBLAS_DISPATCH(cgemm) -> cgemm_incopy
#define	CGEMM_ITCOPY		OPENBLAS_DISPATCH(cgemm) -> cgemm_itcopy

#define	CTRMM_OUNUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_ounucopy
#define	CTRMM_OUTUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_outucopy
#define	CTRMM_OLNUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_olnucopy
#define	CTRMM_OLTUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_oltucopy
#define	CTRSM_OUNUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_ounucopy
#define	CTRSM_OUTUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_outucopy
#define	CTRSM_OLNUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_olnucopy
#define	CTRSM_OLTUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_oltucopy

#define	CTRMM_IUNUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_iunucopy
#define	CTRMM_IUTUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_iutucopy
#define	CTRMM_ILNUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_ilnucopy
#define	CTRMM_ILTUCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_iltucopy
#define	CTRSM_IUNUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_iunucopy
#define	CTRSM_IUTUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_iutucopy
#define	CTRSM_ILNUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_ilnucopy
#define	CTRSM_ILTUCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_iltucopy

#define	CTRMM_OUNNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_ounncopy
#define	CTRMM_OUTNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_outncopy
#define	CTRMM_OLNNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_olnncopy
#define	CTRMM_OLTNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_oltncopy
#define	CTRSM_OUNNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_ounncopy
#define	CTRSM_OUTNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_outncopy
#define	CTRSM_OLNNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_olnncopy
#define	CTRSM_OLTNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_oltncopy

#define	CTRMM_IUNNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_iunncopy
#define	CTRMM_IUTNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_iutncopy
#define	CTRMM_ILNNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_ilnncopy
#define	CTRMM_ILTNCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_iltncopy
#define	CSYMM_INCOPY		    OPENBLAS_DISPATCH(csymm) -> csymm_incopy
#define	CSYMM_ITCOPY		    OPENBLAS_DISPATCH(csymm) -> csymm_itcopy
#define	CTRMM_INCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_incopy
#define	CTRMM_ITCOPY		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_itcopy

#define	CTRSM_IUNNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_iunncopy
#define	CTRSM_IUTNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_iutncopy
#define	CTRSM_ILNNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_ilnncopy
#define	CTRSM_ILTNCOPY		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_iltncopy

#define	SME_CGEMM_KERNEL	OPENBLAS_DISPATCH(cgemm) -> sme_cgemm_kernel
#define	CGEMM_BETA		    OPENBLAS_DISPATCH(cgemm) -> cgemm_beta
#define	CGEMM_KERNEL_N		OPENBLAS_DISPATCH(cgemm) -> cgemm_kernel_n
#define	CGEMM_KERNEL_L		OPENBLAS_DISPATCH(cgemm) -> cgemm_kernel_l
#define	CGEMM_KERNEL_R		OPENBLAS_DISPATCH(cgemm) -> cgemm_kernel_r
#define	CGEMM_KERNEL_B		OPENBLAS_DISPATCH(cgemm) -> cgemm_kernel_b

#define	CSYMM_KERNEL_N		OPENBLAS_DISPATCH(csymm) -> csymm_kernel_n
#define	CSYMM_KERNEL_L		OPENBLAS_DISPATCH(csymm) -> csymm_kernel_l
#define	CSYMM_KERNEL_R		OPENBLAS_DISPATCH(csymm) -> csymm_kernel_r
#define	CSYMM_KERNEL_B		OPENBLAS_DISPATCH(csymm) -> csymm_kernel_b
#define	CTRMM_GEMM_KERNEL_N	OPENBLAS_DISPATCH(ctrmm) -> ctrmm_gemm_kernel_n
#define	CTRMM_GEMM_KERNEL_L	OPENBLAS_DISPATCH(ctrmm) -> ctrmm_gemm_kernel_l
#define	CTRMM_GEMM_KERNEL_R	OPENBLAS_DISPATCH(ctrmm) -> ctrmm_gemm_kernel_r
#define	CTRMM_GEMM_KERNEL_B	OPENBLAS_DISPATCH(ctrmm) -> ctrmm_gemm_kernel_b

#define	CTRMM_KERNEL_LN		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_LN
#define	CTRMM_KERNEL_LT		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_LT
#define	CTRMM_KERNEL_LR		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_LR
#define	CTRMM_KERNEL_LC		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_LC
#define	CTRMM_KERNEL_RN		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_RN
#define	CTRMM_KERNEL_RT		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_RT
#define	CTRMM_KERNEL_RR		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_RR
#define	CTRMM_KERNEL_RC		OPENBLAS_DISPATCH(ctrmm) -> ctrmm_kernel_RC

#define	CTRSM_KERNEL_LN		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_LN
#define	CTRSM_KERNEL_LT		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_LT
#define	CTRSM_KERNEL_LR		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_LR
#define	CTRSM_KERNEL_LC		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_LC
#define	CTRSM_KERNEL_RN		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_RN
#define	CTRSM_KERNEL_RT		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_RT
#define	CTRSM_KERNEL_RR		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_RR
#define	CTRSM_KERNEL_RC		OPENBLAS_DISPATCH(ctrsm) -> ctrsm_kernel_RC

#define	CSYMM_IUTCOPY		OPENBLAS_DISPATCH(csymm) -> csymm_iutcopy
#define	CSYMM_ILTCOPY		OPENBLAS_DISPATCH(csymm) -> csymm_iltcopy
#define	CSYMM_OUTCOPY		OPENBLAS_DISPATCH(csymm) -> csymm_outcopy
#define	CSYMM_OLTCOPY		OPENBLAS_DISPATCH(csymm) -> csymm_oltcopy

#define	CHEMM_OUTCOPY		OPENBLAS_DISPATCH(chemm) -> chemm_outcopy
#define	CHEMM_OLTCOPY		OPENBLAS_DISPATCH(chemm) -> chemm_oltcopy
#define	CHEMM_IUTCOPY		OPENBLAS_DISPATCH(chemm) -> chemm_iutcopy
#define	CHEMM_ILTCOPY		OPENBLAS_DISPATCH(chemm) -> chemm_iltcopy

#define	CGEMM3M_ONCOPYB		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_oncopyb
#define	CGEMM3M_ONCOPYR		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_oncopyr
#define	CGEMM3M_ONCOPYI		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_oncopyi
#define	CGEMM3M_OTCOPYB		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_otcopyb
#define	CGEMM3M_OTCOPYR		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_otcopyr
#define	CGEMM3M_OTCOPYI		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_otcopyi

#define	CGEMM3M_INCOPYB		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_incopyb
#define	CGEMM3M_INCOPYR		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_incopyr
#define	CGEMM3M_INCOPYI		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_incopyi
#define	CGEMM3M_ITCOPYB		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_itcopyb
#define	CGEMM3M_ITCOPYR		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_itcopyr
#define	CGEMM3M_ITCOPYI		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_itcopyi

#define	CSYMM3M_ILCOPYB		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_ilcopyb
#define	CSYMM3M_IUCOPYB		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_iucopyb
#define	CSYMM3M_ILCOPYR		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_ilcopyr
#define	CSYMM3M_IUCOPYR		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_iucopyr
#define	CSYMM3M_ILCOPYI		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_ilcopyi
#define	CSYMM3M_IUCOPYI		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_iucopyi

#define	CSYMM3M_OLCOPYB		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_olcopyb
#define	CSYMM3M_OUCOPYB		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_oucopyb
#define	CSYMM3M_OLCOPYR		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_olcopyr
#define	CSYMM3M_OUCOPYR		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_oucopyr
#define	CSYMM3M_OLCOPYI		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_olcopyi
#define	CSYMM3M_OUCOPYI		OPENBLAS_DISPATCH(csymm3m) -> csymm3m_oucopyi

#define	CHEMM3M_ILCOPYB		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_ilcopyb
#define	CHEMM3M_IUCOPYB		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_iucopyb
#define	CHEMM3M_ILCOPYR		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_ilcopyr
#define	CHEMM3M_IUCOPYR		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_iucopyr
#define	CHEMM3M_ILCOPYI		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_ilcopyi
#define	CHEMM3M_IUCOPYI		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_iucopyi

#define	CHEMM3M_OLCOPYB		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_olcopyb
#define	CHEMM3M_OUCOPYB		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_oucopyb
#define	CHEMM3M_OLCOPYR		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_olcopyr
#define	CHEMM3M_OUCOPYR		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_oucopyr
#define	CHEMM3M_OLCOPYI		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_olcopyi
#define	CHEMM3M_OUCOPYI		OPENBLAS_DISPATCH(chemm3m) -> chemm3m_oucopyi

#define	CGEMM3M_KERNEL		OPENBLAS_DISPATCH(cgemm3m) -> cgemm3m_kernel

#define CNEG_TCOPY		OPENBLAS_DISPATCH(cneg) -> cneg_tcopy
#define CLASWP_NCOPY		OPENBLAS_DISPATCH(claswp) -> claswp_ncopy

#define CAXPBY_K                OPENBLAS_DISPATCH(caxpby) -> caxpby_k

#define COMATCOPY_K_CN          OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_cn
#define COMATCOPY_K_RN          OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_rn
#define COMATCOPY_K_CT          OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_ct
#define COMATCOPY_K_RT          OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_rt
#define COMATCOPY_K_CNC         OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_cnc
#define COMATCOPY_K_RNC         OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_rnc
#define COMATCOPY_K_CTC         OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_ctc
#define COMATCOPY_K_RTC         OPENBLAS_DISPATCH(comatcopy) -> comatcopy_k_rtc

#define CIMATCOPY_K_CN          OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_cn
#define CIMATCOPY_K_RN          OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_rn
#define CIMATCOPY_K_CT          OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_ct
#define CIMATCOPY_K_RT          OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_rt
#define CIMATCOPY_K_CNC         OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_cnc
#define CIMATCOPY_K_RNC         OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_rnc
#define CIMATCOPY_K_CTC         OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_ctc
#define CIMATCOPY_K_RTC         OPENBLAS_DISPATCH(cimatcopy) -> cimatcopy_k_rtc

#define CGEADD_K                OPENBLAS_DISPATCH(cgeadd) -> cgeadd_k 

#define CGEMM_SMALL_MATRIX_PERMIT	OPENBLAS_DISPATCH(cgemm) -> cgemm_small_matrix_permit

#endif

#define CGEMM_SMALL_KERNEL_BASE		OPENBLAS_DISPATCH_BASE(cgemm)
#define CGEMM_SMALL_KERNEL_NN		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_nn)
#define CGEMM_SMALL_KERNEL_NT		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_nt)
#define CGEMM_SMALL_KERNEL_NR		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_nr)
#define CGEMM_SMALL_KERNEL_NC		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_nc)

#define CGEMM_SMALL_KERNEL_TN		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_tn)
#define CGEMM_SMALL_KERNEL_TT		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_tt)
#define CGEMM_SMALL_KERNEL_TR		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_tr)
#define CGEMM_SMALL_KERNEL_TC		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_tc)

#define CGEMM_SMALL_KERNEL_RN		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_rn)
#define CGEMM_SMALL_KERNEL_RT		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_rt)
#define CGEMM_SMALL_KERNEL_RR		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_rr)
#define CGEMM_SMALL_KERNEL_RC		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_rc)

#define CGEMM_SMALL_KERNEL_CN		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_cn)
#define CGEMM_SMALL_KERNEL_CT		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_ct)
#define CGEMM_SMALL_KERNEL_CR		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_cr)
#define CGEMM_SMALL_KERNEL_CC		OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_cc)

#define CGEMM_SMALL_KERNEL_B0_NN	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_nn)
#define CGEMM_SMALL_KERNEL_B0_NT	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_nt)
#define CGEMM_SMALL_KERNEL_B0_NR	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_nr)
#define CGEMM_SMALL_KERNEL_B0_NC	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_nc)

#define CGEMM_SMALL_KERNEL_B0_TN	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_tn)
#define CGEMM_SMALL_KERNEL_B0_TT	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_tt)
#define CGEMM_SMALL_KERNEL_B0_TR	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_tr)
#define CGEMM_SMALL_KERNEL_B0_TC	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_tc)

#define CGEMM_SMALL_KERNEL_B0_RN	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_rn)
#define CGEMM_SMALL_KERNEL_B0_RT	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_rt)
#define CGEMM_SMALL_KERNEL_B0_RR	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_rr)
#define CGEMM_SMALL_KERNEL_B0_RC	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_rc)

#define CGEMM_SMALL_KERNEL_B0_CN	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_cn)
#define CGEMM_SMALL_KERNEL_B0_CT	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_ct)
#define CGEMM_SMALL_KERNEL_B0_CR	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_cr)
#define CGEMM_SMALL_KERNEL_B0_CC	OPENBLAS_DISPATCH_OFFSET(cgemm, cgemm_small_kernel_b0_cc)


#define	CGEMM_NN		cgemm_nn
#define	CGEMM_CN		cgemm_cn
#define	CGEMM_TN		cgemm_tn
#define	CGEMM_NC		cgemm_nc
#define	CGEMM_NT		cgemm_nt
#define	CGEMM_CC		cgemm_cc
#define	CGEMM_CT		cgemm_ct
#define	CGEMM_TC		cgemm_tc
#define	CGEMM_TT		cgemm_tt
#define	CGEMM_NR		cgemm_nr
#define	CGEMM_TR		cgemm_tr
#define	CGEMM_CR		cgemm_cr
#define	CGEMM_RN		cgemm_rn
#define	CGEMM_RT		cgemm_rt
#define	CGEMM_RC		cgemm_rc
#define	CGEMM_RR		cgemm_rr

#define	CSYMM_LU		csymm_LU
#define	CSYMM_LL		csymm_LL
#define	CSYMM_RU		csymm_RU
#define	CSYMM_RL		csymm_RL

#define	CHEMM_LU		chemm_LU
#define	CHEMM_LL		chemm_LL
#define	CHEMM_RU		chemm_RU
#define	CHEMM_RL		chemm_RL

#define	CSYRK_UN		csyrk_UN
#define	CSYRK_UT		csyrk_UT
#define	CSYRK_LN		csyrk_LN
#define	CSYRK_LT		csyrk_LT
#define	CSYRK_UR		csyrk_UN
#define	CSYRK_UC		csyrk_UT
#define	CSYRK_LR		csyrk_LN
#define	CSYRK_LC		csyrk_LT

#define	CSYRK_KERNEL_U		csyrk_kernel_U
#define	CSYRK_KERNEL_L		csyrk_kernel_L

#define	CGEMMT_UNN		cgemmt_UNN
#define	CGEMMT_UNT		cgemmt_UNT
#define	CGEMMT_UNR		cgemmt_UNR
#define	CGEMMT_UNC		cgemmt_UNC
#define	CGEMMT_UTN		cgemmt_UTN
#define	CGEMMT_UTT		cgemmt_UTT
#define	CGEMMT_UTR		cgemmt_UTR
#define	CGEMMT_UTC		cgemmt_UTC
#define	CGEMMT_URN		cgemmt_URN
#define	CGEMMT_URT		cgemmt_URT
#define	CGEMMT_URR		cgemmt_URR
#define	CGEMMT_URC		cgemmt_URC
#define	CGEMMT_UCN		cgemmt_UCN
#define	CGEMMT_UCT		cgemmt_UCT
#define	CGEMMT_UCR		cgemmt_UCR
#define	CGEMMT_UCC		cgemmt_UCC
#define	CGEMMT_LNN		cgemmt_LNN
#define	CGEMMT_LNT		cgemmt_LNT
#define	CGEMMT_LNR		cgemmt_LNR
#define	CGEMMT_LNC		cgemmt_LNC
#define	CGEMMT_LTN		cgemmt_LTN
#define	CGEMMT_LTT		cgemmt_LTT
#define	CGEMMT_LTR		cgemmt_LTR
#define	CGEMMT_LTC		cgemmt_LTC
#define	CGEMMT_LRN		cgemmt_LRN
#define	CGEMMT_LRT		cgemmt_LRT
#define	CGEMMT_LRR		cgemmt_LRR
#define	CGEMMT_LRC		cgemmt_LRC
#define	CGEMMT_LCN		cgemmt_LCN
#define	CGEMMT_LCT		cgemmt_LCT
#define	CGEMMT_LCR		cgemmt_LCR
#define	CGEMMT_LCC		cgemmt_LCC

#define	CGEMMT_KERNEL_UCN	cgemmt_kernel_UCN
#define	CGEMMT_KERNEL_UNC	cgemmt_kernel_UNC
#define	CGEMMT_KERNEL_UCC	cgemmt_kernel_UCC
#define	CGEMMT_KERNEL_LCN	cgemmt_kernel_LCN
#define	CGEMMT_KERNEL_LNC	cgemmt_kernel_LNC
#define	CGEMMT_KERNEL_LCC	cgemmt_kernel_LCC

#define	CHERK_UN		cherk_UN
#define	CHERK_LN		cherk_LN
#define	CHERK_UC		cherk_UC
#define	CHERK_LC		cherk_LC

#define	CHER2K_UN		cher2k_UN
#define	CHER2K_LN		cher2k_LN
#define	CHER2K_UC		cher2k_UC
#define	CHER2K_LC		cher2k_LC

#define	CSYR2K_UN		csyr2k_UN
#define	CSYR2K_UT		csyr2k_UT
#define	CSYR2K_LN		csyr2k_LN
#define	CSYR2K_LT		csyr2k_LT
#define	CSYR2K_UR		csyr2k_UN
#define	CSYR2K_UC		csyr2k_UT
#define	CSYR2K_LR		csyr2k_LN
#define	CSYR2K_LC		csyr2k_LT

#define	CSYR2K_KERNEL_U		csyr2k_kernel_U
#define	CSYR2K_KERNEL_L		csyr2k_kernel_L

#define	CTRMM_LNUU		ctrmm_LNUU
#define	CTRMM_LNUN		ctrmm_LNUN
#define	CTRMM_LNLU		ctrmm_LNLU
#define	CTRMM_LNLN		ctrmm_LNLN
#define	CTRMM_LTUU		ctrmm_LTUU
#define	CTRMM_LTUN		ctrmm_LTUN
#define	CTRMM_LTLU		ctrmm_LTLU
#define	CTRMM_LTLN		ctrmm_LTLN
#define	CTRMM_LRUU		ctrmm_LRUU
#define	CTRMM_LRUN		ctrmm_LRUN
#define	CTRMM_LRLU		ctrmm_LRLU
#define	CTRMM_LRLN		ctrmm_LRLN
#define	CTRMM_LCUU		ctrmm_LCUU
#define	CTRMM_LCUN		ctrmm_LCUN
#define	CTRMM_LCLU		ctrmm_LCLU
#define	CTRMM_LCLN		ctrmm_LCLN
#define	CTRMM_RNUU		ctrmm_RNUU
#define	CTRMM_RNUN		ctrmm_RNUN
#define	CTRMM_RNLU		ctrmm_RNLU
#define	CTRMM_RNLN		ctrmm_RNLN
#define	CTRMM_RTUU		ctrmm_RTUU
#define	CTRMM_RTUN		ctrmm_RTUN
#define	CTRMM_RTLU		ctrmm_RTLU
#define	CTRMM_RTLN		ctrmm_RTLN
#define	CTRMM_RRUU		ctrmm_RRUU
#define	CTRMM_RRUN		ctrmm_RRUN
#define	CTRMM_RRLU		ctrmm_RRLU
#define	CTRMM_RRLN		ctrmm_RRLN
#define	CTRMM_RCUU		ctrmm_RCUU
#define	CTRMM_RCUN		ctrmm_RCUN
#define	CTRMM_RCLU		ctrmm_RCLU
#define	CTRMM_RCLN		ctrmm_RCLN

#define	CTRSM_LNUU		ctrsm_LNUU
#define	CTRSM_LNUN		ctrsm_LNUN
#define	CTRSM_LNLU		ctrsm_LNLU
#define	CTRSM_LNLN		ctrsm_LNLN
#define	CTRSM_LTUU		ctrsm_LTUU
#define	CTRSM_LTUN		ctrsm_LTUN
#define	CTRSM_LTLU		ctrsm_LTLU
#define	CTRSM_LTLN		ctrsm_LTLN
#define	CTRSM_LRUU		ctrsm_LRUU
#define	CTRSM_LRUN		ctrsm_LRUN
#define	CTRSM_LRLU		ctrsm_LRLU
#define	CTRSM_LRLN		ctrsm_LRLN
#define	CTRSM_LCUU		ctrsm_LCUU
#define	CTRSM_LCUN		ctrsm_LCUN
#define	CTRSM_LCLU		ctrsm_LCLU
#define	CTRSM_LCLN		ctrsm_LCLN
#define	CTRSM_RNUU		ctrsm_RNUU
#define	CTRSM_RNUN		ctrsm_RNUN
#define	CTRSM_RNLU		ctrsm_RNLU
#define	CTRSM_RNLN		ctrsm_RNLN
#define	CTRSM_RTUU		ctrsm_RTUU
#define	CTRSM_RTUN		ctrsm_RTUN
#define	CTRSM_RTLU		ctrsm_RTLU
#define	CTRSM_RTLN		ctrsm_RTLN
#define	CTRSM_RRUU		ctrsm_RRUU
#define	CTRSM_RRUN		ctrsm_RRUN
#define	CTRSM_RRLU		ctrsm_RRLU
#define	CTRSM_RRLN		ctrsm_RRLN
#define	CTRSM_RCUU		ctrsm_RCUU
#define	CTRSM_RCUN		ctrsm_RCUN
#define	CTRSM_RCLU		ctrsm_RCLU
#define	CTRSM_RCLN		ctrsm_RCLN

#define	CGEMM_THREAD_NN		cgemm_thread_nn
#define	CGEMM_THREAD_CN		cgemm_thread_cn
#define	CGEMM_THREAD_TN		cgemm_thread_tn
#define	CGEMM_THREAD_NC		cgemm_thread_nc
#define	CGEMM_THREAD_NT		cgemm_thread_nt
#define	CGEMM_THREAD_CC		cgemm_thread_cc
#define	CGEMM_THREAD_CT		cgemm_thread_ct
#define	CGEMM_THREAD_TC		cgemm_thread_tc
#define	CGEMM_THREAD_TT		cgemm_thread_tt
#define	CGEMM_THREAD_NR		cgemm_thread_nr
#define	CGEMM_THREAD_TR		cgemm_thread_tr
#define	CGEMM_THREAD_CR		cgemm_thread_cr
#define	CGEMM_THREAD_RN		cgemm_thread_rn
#define	CGEMM_THREAD_RT		cgemm_thread_rt
#define	CGEMM_THREAD_RC		cgemm_thread_rc
#define	CGEMM_THREAD_RR		cgemm_thread_rr

#define	CSYMM_THREAD_LU		csymm_thread_LU
#define	CSYMM_THREAD_LL		csymm_thread_LL
#define	CSYMM_THREAD_RU		csymm_thread_RU
#define	CSYMM_THREAD_RL		csymm_thread_RL

#define	CHEMM_THREAD_LU		chemm_thread_LU
#define	CHEMM_THREAD_LL		chemm_thread_LL
#define	CHEMM_THREAD_RU		chemm_thread_RU
#define	CHEMM_THREAD_RL		chemm_thread_RL

#define	CSYRK_THREAD_UN		csyrk_thread_UN
#define	CSYRK_THREAD_UT		csyrk_thread_UT
#define	CSYRK_THREAD_LN		csyrk_thread_LN
#define	CSYRK_THREAD_LT		csyrk_thread_LT
#define	CSYRK_THREAD_UR		csyrk_thread_UN
#define	CSYRK_THREAD_UC		csyrk_thread_UT
#define	CSYRK_THREAD_LR		csyrk_thread_LN
#define	CSYRK_THREAD_LC		csyrk_thread_LT

#define	CHERK_THREAD_UN		cherk_thread_UN
#define	CHERK_THREAD_UT		cherk_thread_UT
#define	CHERK_THREAD_LN		cherk_thread_LN
#define	CHERK_THREAD_LT		cherk_thread_LT
#define	CHERK_THREAD_UR		cherk_thread_UR
#define	CHERK_THREAD_UC		cherk_thread_UC
#define	CHERK_THREAD_LR		cherk_thread_LR
#define	CHERK_THREAD_LC		cherk_thread_LC

#define	CGEMM3M_NN		cgemm3m_nn
#define	CGEMM3M_CN		cgemm3m_cn
#define	CGEMM3M_TN		cgemm3m_tn
#define	CGEMM3M_NC		cgemm3m_nc
#define	CGEMM3M_NT		cgemm3m_nt
#define	CGEMM3M_CC		cgemm3m_cc
#define	CGEMM3M_CT		cgemm3m_ct
#define	CGEMM3M_TC		cgemm3m_tc
#define	CGEMM3M_TT		cgemm3m_tt
#define	CGEMM3M_NR		cgemm3m_nr
#define	CGEMM3M_TR		cgemm3m_tr
#define	CGEMM3M_CR		cgemm3m_cr
#define	CGEMM3M_RN		cgemm3m_rn
#define	CGEMM3M_RT		cgemm3m_rt
#define	CGEMM3M_RC		cgemm3m_rc
#define	CGEMM3M_RR		cgemm3m_rr

#define	CGEMM3M_THREAD_NN	cgemm3m_thread_nn
#define	CGEMM3M_THREAD_CN	cgemm3m_thread_cn
#define	CGEMM3M_THREAD_TN	cgemm3m_thread_tn
#define	CGEMM3M_THREAD_NC	cgemm3m_thread_nc
#define	CGEMM3M_THREAD_NT	cgemm3m_thread_nt
#define	CGEMM3M_THREAD_CC	cgemm3m_thread_cc
#define	CGEMM3M_THREAD_CT	cgemm3m_thread_ct
#define	CGEMM3M_THREAD_TC	cgemm3m_thread_tc
#define	CGEMM3M_THREAD_TT	cgemm3m_thread_tt
#define	CGEMM3M_THREAD_NR	cgemm3m_thread_nr
#define	CGEMM3M_THREAD_TR	cgemm3m_thread_tr
#define	CGEMM3M_THREAD_CR	cgemm3m_thread_cr
#define	CGEMM3M_THREAD_RN	cgemm3m_thread_rn
#define	CGEMM3M_THREAD_RT	cgemm3m_thread_rt
#define	CGEMM3M_THREAD_RC	cgemm3m_thread_rc
#define	CGEMM3M_THREAD_RR	cgemm3m_thread_rr

#define	CSYMM3M_LU		csymm3m_LU
#define	CSYMM3M_LL		csymm3m_LL
#define	CSYMM3M_RU		csymm3m_RU
#define	CSYMM3M_RL		csymm3m_RL

#define	CSYMM3M_THREAD_LU	csymm3m_thread_LU
#define	CSYMM3M_THREAD_LL	csymm3m_thread_LL
#define	CSYMM3M_THREAD_RU	csymm3m_thread_RU
#define	CSYMM3M_THREAD_RL	csymm3m_thread_RL

#define	CHEMM3M_LU		chemm3m_LU
#define	CHEMM3M_LL		chemm3m_LL
#define	CHEMM3M_RU		chemm3m_RU
#define	CHEMM3M_RL		chemm3m_RL

#define	CHEMM3M_THREAD_LU	chemm3m_thread_LU
#define	CHEMM3M_THREAD_LL	chemm3m_thread_LL
#define	CHEMM3M_THREAD_RU	chemm3m_thread_RU
#define	CHEMM3M_THREAD_RL	chemm3m_thread_RL

#endif
