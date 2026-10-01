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

#ifndef COMMON_D_H
#define COMMON_D_H

#ifndef DYNAMIC_ARCH

#define	DAMAX_K			damax_k
#define	DAMIN_K			damin_k
#define	DMAX_K			dmax_k
#define	DMIN_K			dmin_k
#define	IDAMAX_K		idamax_k
#define	IDAMIN_K		idamin_k
#define	IDMAX_K			idmax_k
#define	IDMIN_K			idmin_k
#define	DASUM_K			dasum_k
#define	DAXPYU_K		daxpy_k
#define	DAXPYC_K		daxpy_k
#define	DCOPY_K			dcopy_k
#define	DDOTU_K			ddot_k
#define	DDOTC_K			ddot_k
#define	DNRM2_K			dnrm2_k
#define	DSCAL_K			dscal_k
#define	DSUM_K			dsum_k
#define	DSWAP_K			dswap_k
#define	DROT_K			drot_k
#define DROTM_K         drotm_k

#define	DGEMV_N			dgemv_n
#define	DGEMV_T			dgemv_t
#define	DGEMV_R			dgemv_n
#define	DGEMV_C			dgemv_t
#define	DGEMV_O			dgemv_n
#define	DGEMV_U			dgemv_t
#define	DGEMV_S			dgemv_n
#define	DGEMV_D			dgemv_t

#define	DGERU_K			dger_k
#define	DGERC_K			dger_k
#define	DGERV_K			dger_k
#define	DGERD_K			dger_k

#define DSYMV_U			dsymv_U
#define DSYMV_L			dsymv_L

#define DSYMV_THREAD_U		dsymv_thread_U
#define DSYMV_THREAD_L		dsymv_thread_L

#define	DGEMM_ONCOPY		dgemm_oncopy
#define	DGEMM_OTCOPY		dgemm_otcopy

#if DGEMM_DEFAULT_UNROLL_M == DGEMM_DEFAULT_UNROLL_N
#define	DGEMM_INCOPY		dgemm_oncopy
#define	DGEMM_ITCOPY		dgemm_otcopy
#else
#define	DGEMM_INCOPY		dgemm_incopy
#define	DGEMM_ITCOPY		dgemm_itcopy
#endif

#define	DSYMM_INCOPY	    	dsymm_incopy
#define	DSYMM_ITCOPY	    	dsymm_itcopy
#define	DTRMM_INCOPY		dtrmm_incopy
#define	DTRMM_ITCOPY		dtrmm_itcopy

#define	DTRMM_OUNUCOPY		dtrmm_ounucopy
#define	DTRMM_OUNNCOPY		dtrmm_ounncopy
#define	DTRMM_OUTUCOPY		dtrmm_outucopy
#define	DTRMM_OUTNCOPY		dtrmm_outncopy
#define	DTRMM_OLNUCOPY		dtrmm_olnucopy
#define	DTRMM_OLNNCOPY		dtrmm_olnncopy
#define	DTRMM_OLTUCOPY		dtrmm_oltucopy
#define	DTRMM_OLTNCOPY		dtrmm_oltncopy

#define	DTRSM_OUNUCOPY		dtrsm_ounucopy
#define	DTRSM_OUNNCOPY		dtrsm_ounncopy
#define	DTRSM_OUTUCOPY		dtrsm_outucopy
#define	DTRSM_OUTNCOPY		dtrsm_outncopy
#define	DTRSM_OLNUCOPY		dtrsm_olnucopy
#define	DTRSM_OLNNCOPY		dtrsm_olnncopy
#define	DTRSM_OLTUCOPY		dtrsm_oltucopy
#define	DTRSM_OLTNCOPY		dtrsm_oltncopy

#if DGEMM_DEFAULT_UNROLL_M == DGEMM_DEFAULT_UNROLL_N
#define	DTRMM_IUNUCOPY		dtrmm_ounucopy
#define	DTRMM_IUNNCOPY		dtrmm_ounncopy
#define	DTRMM_IUTUCOPY		dtrmm_outucopy
#define	DTRMM_IUTNCOPY		dtrmm_outncopy
#define	DTRMM_ILNUCOPY		dtrmm_olnucopy
#define	DTRMM_ILNNCOPY		dtrmm_olnncopy
#define	DTRMM_ILTUCOPY		dtrmm_oltucopy
#define	DTRMM_ILTNCOPY		dtrmm_oltncopy

#define	DTRSM_IUNUCOPY		dtrsm_ounucopy
#define	DTRSM_IUNNCOPY		dtrsm_ounncopy
#define	DTRSM_IUTUCOPY		dtrsm_outucopy
#define	DTRSM_IUTNCOPY		dtrsm_outncopy
#define	DTRSM_ILNUCOPY		dtrsm_olnucopy
#define	DTRSM_ILNNCOPY		dtrsm_olnncopy
#define	DTRSM_ILTUCOPY		dtrsm_oltucopy
#define	DTRSM_ILTNCOPY		dtrsm_oltncopy
#else
#define	DTRMM_IUNUCOPY		dtrmm_iunucopy
#define	DTRMM_IUNNCOPY		dtrmm_iunncopy
#define	DTRMM_IUTUCOPY		dtrmm_iutucopy
#define	DTRMM_IUTNCOPY		dtrmm_iutncopy
#define	DTRMM_ILNUCOPY		dtrmm_ilnucopy
#define	DTRMM_ILNNCOPY		dtrmm_ilnncopy
#define	DTRMM_ILTUCOPY		dtrmm_iltucopy
#define	DTRMM_ILTNCOPY		dtrmm_iltncopy

#define	DTRSM_IUNUCOPY		dtrsm_iunucopy
#define	DTRSM_IUNNCOPY		dtrsm_iunncopy
#define	DTRSM_IUTUCOPY		dtrsm_iutucopy
#define	DTRSM_IUTNCOPY		dtrsm_iutncopy
#define	DTRSM_ILNUCOPY		dtrsm_ilnucopy
#define	DTRSM_ILNNCOPY		dtrsm_ilnncopy
#define	DTRSM_ILTUCOPY		dtrsm_iltucopy
#define	DTRSM_ILTNCOPY		dtrsm_iltncopy
#endif

#define	DGEMM_BETA		dgemm_beta

#define	DGEMM_KERNEL		dgemm_kernel
#define SME_DGEMM_KERNEL	sme_dgemm_kernel
#define	DSYMM_KERNEL		dsymm_kernel
#define	DTRMM_GEMM_KERNEL	dtrmm_gemm_kernel

#define	DTRMM_KERNEL_LN		dtrmm_kernel_LN
#define	DTRMM_KERNEL_LT		dtrmm_kernel_LT
#define	DTRMM_KERNEL_LR		dtrmm_kernel_LN
#define	DTRMM_KERNEL_LC		dtrmm_kernel_LT
#define	DTRMM_KERNEL_RN		dtrmm_kernel_RN
#define	DTRMM_KERNEL_RT		dtrmm_kernel_RT
#define	DTRMM_KERNEL_RR		dtrmm_kernel_RN
#define	DTRMM_KERNEL_RC		dtrmm_kernel_RT

#define	DTRSM_KERNEL_LN		dtrsm_kernel_LN
#define	DTRSM_KERNEL_LT		dtrsm_kernel_LT
#define	DTRSM_KERNEL_LR		dtrsm_kernel_LN
#define	DTRSM_KERNEL_LC		dtrsm_kernel_LT
#define	DTRSM_KERNEL_RN		dtrsm_kernel_RN
#define	DTRSM_KERNEL_RT		dtrsm_kernel_RT
#define	DTRSM_KERNEL_RR		dtrsm_kernel_RN
#define	DTRSM_KERNEL_RC		dtrsm_kernel_RT

#define	DSYMM_OUTCOPY		dsymm_outcopy
#define	DSYMM_OLTCOPY		dsymm_oltcopy
#if DGEMM_DEFAULT_UNROLL_M == DGEMM_DEFAULT_UNROLL_N
#define	DSYMM_IUTCOPY		dsymm_outcopy
#define	DSYMM_ILTCOPY		dsymm_oltcopy
#else
#define	DSYMM_IUTCOPY		dsymm_iutcopy
#define	DSYMM_ILTCOPY		dsymm_iltcopy
#endif

#define DNEG_TCOPY		dneg_tcopy
#define DLASWP_NCOPY		dlaswp_ncopy

#define	DAXPBY_K		daxpby_k
#define DOMATCOPY_K_CN		domatcopy_k_cn
#define DOMATCOPY_K_RN		domatcopy_k_rn
#define DOMATCOPY_K_CT		domatcopy_k_ct
#define DOMATCOPY_K_RT		domatcopy_k_rt

#define DIMATCOPY_K_CN		dimatcopy_k_cn
#define DIMATCOPY_K_RN		dimatcopy_k_rn
#define DIMATCOPY_K_CT      dimatcopy_k_ct
#define DIMATCOPY_K_RT      dimatcopy_k_rt
#define DGEADD_K                dgeadd_k 

#define DGEMM_SMALL_MATRIX_PERMIT	dgemm_small_matrix_permit

#else

#define	DAMAX_K			OPENBLAS_DISPATCH(damax) -> damax_k
#define	DAMIN_K			OPENBLAS_DISPATCH(damin) -> damin_k
#define	DMAX_K			OPENBLAS_DISPATCH(dmax) -> dmax_k
#define	DMIN_K			OPENBLAS_DISPATCH(dmin) -> dmin_k
#define	IDAMAX_K		OPENBLAS_DISPATCH(idamax) -> idamax_k
#define	IDAMIN_K		OPENBLAS_DISPATCH(idamin) -> idamin_k
#define	IDMAX_K			OPENBLAS_DISPATCH(idmax) -> idmax_k
#define	IDMIN_K			OPENBLAS_DISPATCH(idmin) -> idmin_k
#define	DASUM_K			OPENBLAS_DISPATCH(dasum) -> dasum_k
#define	DAXPYU_K		OPENBLAS_DISPATCH(daxpy) -> daxpy_k
#define	DAXPYC_K		OPENBLAS_DISPATCH(daxpy) -> daxpy_k
#define	DCOPY_K			OPENBLAS_DISPATCH(dcopy) -> dcopy_k
#define	DDOTU_K			OPENBLAS_DISPATCH(ddot) -> ddot_k
#define	DDOTC_K			OPENBLAS_DISPATCH(ddot) -> ddot_k
#define	DNRM2_K			OPENBLAS_DISPATCH(dnrm2) -> dnrm2_k
#define	DSCAL_K			OPENBLAS_DISPATCH(dscal) -> dscal_k
#define	DSUM_K			OPENBLAS_DISPATCH(dsum) -> dsum_k
#define	DSWAP_K			OPENBLAS_DISPATCH(dswap) -> dswap_k
#define	DROT_K			OPENBLAS_DISPATCH(drot) -> drot_k
#define	DROTM_K			OPENBLAS_DISPATCH(drotm) -> drotm_k

#define	DGEMV_N			OPENBLAS_DISPATCH(dgemv) -> dgemv_n
#define	DGEMV_T			OPENBLAS_DISPATCH(dgemv) -> dgemv_t
#define	DGEMV_R			OPENBLAS_DISPATCH(dgemv) -> dgemv_n
#define	DGEMV_C			OPENBLAS_DISPATCH(dgemv) -> dgemv_t
#define	DGEMV_O			OPENBLAS_DISPATCH(dgemv) -> dgemv_n
#define	DGEMV_U			OPENBLAS_DISPATCH(dgemv) -> dgemv_t
#define	DGEMV_S			OPENBLAS_DISPATCH(dgemv) -> dgemv_n
#define	DGEMV_D			OPENBLAS_DISPATCH(dgemv) -> dgemv_t

#define	DGERU_K			OPENBLAS_DISPATCH(dger) -> dger_k
#define	DGERC_K			OPENBLAS_DISPATCH(dger) -> dger_k
#define	DGERV_K			OPENBLAS_DISPATCH(dger) -> dger_k
#define	DGERD_K			OPENBLAS_DISPATCH(dger) -> dger_k

#define DSYMV_U			OPENBLAS_DISPATCH(dsymv) -> dsymv_U
#define DSYMV_L			OPENBLAS_DISPATCH(dsymv) -> dsymv_L

#define DSYMV_THREAD_U		dsymv_thread_U
#define DSYMV_THREAD_L		dsymv_thread_L

#define	DGEMM_ONCOPY		OPENBLAS_DISPATCH(dgemm) -> dgemm_oncopy
#define	DGEMM_OTCOPY		OPENBLAS_DISPATCH(dgemm) -> dgemm_otcopy
#define	DGEMM_INCOPY		OPENBLAS_DISPATCH(dgemm) -> dgemm_incopy
#define	DGEMM_ITCOPY		OPENBLAS_DISPATCH(dgemm) -> dgemm_itcopy

#define	DTRMM_OUNUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_ounucopy
#define	DTRMM_OUTUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_outucopy
#define	DTRMM_OLNUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_olnucopy
#define	DTRMM_OLTUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_oltucopy
#define	DTRSM_OUNUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_ounucopy
#define	DTRSM_OUTUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_outucopy
#define	DTRSM_OLNUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_olnucopy
#define	DTRSM_OLTUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_oltucopy

#define	DTRMM_IUNUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_iunucopy
#define	DTRMM_IUTUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_iutucopy
#define	DTRMM_ILNUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_ilnucopy
#define	DTRMM_ILTUCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_iltucopy
#define	DTRSM_IUNUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_iunucopy
#define	DTRSM_IUTUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_iutucopy
#define	DTRSM_ILNUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_ilnucopy
#define	DTRSM_ILTUCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_iltucopy

#define	DTRMM_OUNNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_ounncopy
#define	DTRMM_OUTNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_outncopy
#define	DTRMM_OLNNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_olnncopy
#define	DTRMM_OLTNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_oltncopy
#define	DTRSM_OUNNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_ounncopy
#define	DTRSM_OUTNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_outncopy
#define	DTRSM_OLNNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_olnncopy
#define	DTRSM_OLTNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_oltncopy

#define	DTRMM_IUNNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_iunncopy
#define	DTRMM_IUTNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_iutncopy
#define	DTRMM_ILNNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_ilnncopy
#define	DTRMM_ILTNCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_iltncopy
#define	DSYMM_INCOPY		    OPENBLAS_DISPATCH(dsymm) -> dsymm_incopy
#define	DSYMM_ITCOPY		    OPENBLAS_DISPATCH(dsymm) -> dsymm_itcopy
#define	DTRMM_INCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_incopy
#define	DTRMM_ITCOPY		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_itcopy

#define	DTRSM_IUNNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_iunncopy
#define	DTRSM_IUTNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_iutncopy
#define	DTRSM_ILNNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_ilnncopy
#define	DTRSM_ILTNCOPY		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_iltncopy

#define	DGEMM_BETA	    	OPENBLAS_DISPATCH(dgemm) -> dgemm_beta
#define	DGEMM_KERNEL		OPENBLAS_DISPATCH(dgemm) -> dgemm_kernel
#define	SME_DGEMM_KERNEL	OPENBLAS_DISPATCH(dgemm) -> sme_dgemm_kernel
#define	DSYMM_KERNEL		OPENBLAS_DISPATCH(dsymm) -> dsymm_kernel
#define	DTRMM_GEMM_KERNEL	OPENBLAS_DISPATCH(dtrmm) -> dtrmm_gemm_kernel

#define	DTRMM_KERNEL_LN		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_LN
#define	DTRMM_KERNEL_LT		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_LT
#define	DTRMM_KERNEL_LR		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_LN
#define	DTRMM_KERNEL_LC		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_LT
#define	DTRMM_KERNEL_RN		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_RN
#define	DTRMM_KERNEL_RT		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_RT
#define	DTRMM_KERNEL_RR		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_RN
#define	DTRMM_KERNEL_RC		OPENBLAS_DISPATCH(dtrmm) -> dtrmm_kernel_RT

#define	DTRSM_KERNEL_LN		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_LN
#define	DTRSM_KERNEL_LT		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_LT
#define	DTRSM_KERNEL_LR		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_LN
#define	DTRSM_KERNEL_LC		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_LT
#define	DTRSM_KERNEL_RN		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_RN
#define	DTRSM_KERNEL_RT		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_RT
#define	DTRSM_KERNEL_RR		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_RN
#define	DTRSM_KERNEL_RC		OPENBLAS_DISPATCH(dtrsm) -> dtrsm_kernel_RT

#define	DSYMM_IUTCOPY		OPENBLAS_DISPATCH(dsymm) -> dsymm_iutcopy
#define	DSYMM_ILTCOPY		OPENBLAS_DISPATCH(dsymm) -> dsymm_iltcopy
#define	DSYMM_OUTCOPY		OPENBLAS_DISPATCH(dsymm) -> dsymm_outcopy
#define	DSYMM_OLTCOPY		OPENBLAS_DISPATCH(dsymm) -> dsymm_oltcopy

#define DNEG_TCOPY		OPENBLAS_DISPATCH(dneg) -> dneg_tcopy
#define DLASWP_NCOPY		OPENBLAS_DISPATCH(dlaswp) -> dlaswp_ncopy

#define	DAXPBY_K		OPENBLAS_DISPATCH(daxpby) -> daxpby_k
#define DOMATCOPY_K_CN		OPENBLAS_DISPATCH(domatcopy) -> domatcopy_k_cn
#define DOMATCOPY_K_RN		OPENBLAS_DISPATCH(domatcopy) -> domatcopy_k_rn
#define DOMATCOPY_K_CT		OPENBLAS_DISPATCH(domatcopy) -> domatcopy_k_ct
#define DOMATCOPY_K_RT		OPENBLAS_DISPATCH(domatcopy) -> domatcopy_k_rt
#define DIMATCOPY_K_CN		OPENBLAS_DISPATCH(dimatcopy) -> dimatcopy_k_cn
#define DIMATCOPY_K_RN		OPENBLAS_DISPATCH(dimatcopy) -> dimatcopy_k_rn
#define DIMATCOPY_K_CT		OPENBLAS_DISPATCH(dimatcopy) -> dimatcopy_k_ct
#define DIMATCOPY_K_RT		OPENBLAS_DISPATCH(dimatcopy) -> dimatcopy_k_rt

#define DGEADD_K                OPENBLAS_DISPATCH(dgeadd) -> dgeadd_k 

#define DGEMM_SMALL_MATRIX_PERMIT	OPENBLAS_DISPATCH(dgemm) -> dgemm_small_matrix_permit

#endif

#define DGEMM_SMALL_KERNEL_BASE		OPENBLAS_DISPATCH_BASE(dgemm)
#define DGEMM_SMALL_KERNEL_NN		OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_nn)
#define DGEMM_SMALL_KERNEL_NT		OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_nt)
#define DGEMM_SMALL_KERNEL_TN		OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_tn)
#define DGEMM_SMALL_KERNEL_TT		OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_tt)

#define DGEMM_SMALL_KERNEL_B0_NN	OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_b0_nn)
#define DGEMM_SMALL_KERNEL_B0_NT	OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_b0_nt)
#define DGEMM_SMALL_KERNEL_B0_TN	OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_b0_tn)
#define DGEMM_SMALL_KERNEL_B0_TT	OPENBLAS_DISPATCH_OFFSET(dgemm, dgemm_small_kernel_b0_tt)


#define	DGEMM_NN		dgemm_nn
#define	DGEMM_CN		dgemm_tn
#define	DGEMM_TN		dgemm_tn
#define	DGEMM_NC		dgemm_nt
#define	DGEMM_NT		dgemm_nt
#define	DGEMM_CC		dgemm_tt
#define	DGEMM_CT		dgemm_tt
#define	DGEMM_TC		dgemm_tt
#define	DGEMM_TT		dgemm_tt
#define	DGEMM_NR		dgemm_nn
#define	DGEMM_TR		dgemm_tn
#define	DGEMM_CR		dgemm_tn
#define	DGEMM_RN		dgemm_nn
#define	DGEMM_RT		dgemm_nt
#define	DGEMM_RC		dgemm_nt
#define	DGEMM_RR		dgemm_nn

#define	DSYMM_LU		dsymm_LU
#define	DSYMM_LL		dsymm_LL
#define	DSYMM_RU		dsymm_RU
#define	DSYMM_RL		dsymm_RL

#define	DHEMM_LU		dhemm_LU
#define	DHEMM_LL		dhemm_LL
#define	DHEMM_RU		dhemm_RU
#define	DHEMM_RL		dhemm_RL

#define	DSYRK_UN		dsyrk_UN
#define	DSYRK_UT		dsyrk_UT
#define	DSYRK_LN		dsyrk_LN
#define	DSYRK_LT		dsyrk_LT
#define	DSYRK_UR		dsyrk_UN
#define	DSYRK_UC		dsyrk_UT
#define	DSYRK_LR		dsyrk_LN
#define	DSYRK_LC		dsyrk_LT

#define	DSYRK_KERNEL_U		dsyrk_kernel_U
#define	DSYRK_KERNEL_L		dsyrk_kernel_L

#define	DGEMMT_UNN		dgemmt_UNN
#define	DGEMMT_UNT		dgemmt_UNT
#define	DGEMMT_UTN		dgemmt_UTN
#define	DGEMMT_UTT		dgemmt_UTT
#define	DGEMMT_LNN		dgemmt_LNN
#define	DGEMMT_LNT		dgemmt_LNT
#define	DGEMMT_LTN		dgemmt_LTN
#define	DGEMMT_LTT		dgemmt_LTT

#define	DHERK_UN		dsyrk_UN
#define	DHERK_LN		dsyrk_LN
#define	DHERK_UC		dsyrk_UT
#define	DHERK_LC		dsyrk_LT

#define	DHER2K_UN		dsyr2k_UN
#define	DHER2K_LN		dsyr2k_LN
#define	DHER2K_UC		dsyr2k_UT
#define	DHER2K_LC		dsyr2k_LT

#define	DSYR2K_UN		dsyr2k_UN
#define	DSYR2K_UT		dsyr2k_UT
#define	DSYR2K_LN		dsyr2k_LN
#define	DSYR2K_LT		dsyr2k_LT
#define	DSYR2K_UR		dsyr2k_UN
#define	DSYR2K_UC		dsyr2k_UT
#define	DSYR2K_LR		dsyr2k_LN
#define	DSYR2K_LC		dsyr2k_LT

#define	DSYR2K_KERNEL_U		dsyr2k_kernel_U
#define	DSYR2K_KERNEL_L		dsyr2k_kernel_L

#define	DTRMM_LNUU		dtrmm_LNUU
#define	DTRMM_LNUN		dtrmm_LNUN
#define	DTRMM_LNLU		dtrmm_LNLU
#define	DTRMM_LNLN		dtrmm_LNLN
#define	DTRMM_LTUU		dtrmm_LTUU
#define	DTRMM_LTUN		dtrmm_LTUN
#define	DTRMM_LTLU		dtrmm_LTLU
#define	DTRMM_LTLN		dtrmm_LTLN
#define	DTRMM_LRUU		dtrmm_LNUU
#define	DTRMM_LRUN		dtrmm_LNUN
#define	DTRMM_LRLU		dtrmm_LNLU
#define	DTRMM_LRLN		dtrmm_LNLN
#define	DTRMM_LCUU		dtrmm_LTUU
#define	DTRMM_LCUN		dtrmm_LTUN
#define	DTRMM_LCLU		dtrmm_LTLU
#define	DTRMM_LCLN		dtrmm_LTLN
#define	DTRMM_RNUU		dtrmm_RNUU
#define	DTRMM_RNUN		dtrmm_RNUN
#define	DTRMM_RNLU		dtrmm_RNLU
#define	DTRMM_RNLN		dtrmm_RNLN
#define	DTRMM_RTUU		dtrmm_RTUU
#define	DTRMM_RTUN		dtrmm_RTUN
#define	DTRMM_RTLU		dtrmm_RTLU
#define	DTRMM_RTLN		dtrmm_RTLN
#define	DTRMM_RRUU		dtrmm_RNUU
#define	DTRMM_RRUN		dtrmm_RNUN
#define	DTRMM_RRLU		dtrmm_RNLU
#define	DTRMM_RRLN		dtrmm_RNLN
#define	DTRMM_RCUU		dtrmm_RTUU
#define	DTRMM_RCUN		dtrmm_RTUN
#define	DTRMM_RCLU		dtrmm_RTLU
#define	DTRMM_RCLN		dtrmm_RTLN

#define	DTRSM_LNUU		dtrsm_LNUU
#define	DTRSM_LNUN		dtrsm_LNUN
#define	DTRSM_LNLU		dtrsm_LNLU
#define	DTRSM_LNLN		dtrsm_LNLN
#define	DTRSM_LTUU		dtrsm_LTUU
#define	DTRSM_LTUN		dtrsm_LTUN
#define	DTRSM_LTLU		dtrsm_LTLU
#define	DTRSM_LTLN		dtrsm_LTLN
#define	DTRSM_LRUU		dtrsm_LNUU
#define	DTRSM_LRUN		dtrsm_LNUN
#define	DTRSM_LRLU		dtrsm_LNLU
#define	DTRSM_LRLN		dtrsm_LNLN
#define	DTRSM_LCUU		dtrsm_LTUU
#define	DTRSM_LCUN		dtrsm_LTUN
#define	DTRSM_LCLU		dtrsm_LTLU
#define	DTRSM_LCLN		dtrsm_LTLN
#define	DTRSM_RNUU		dtrsm_RNUU
#define	DTRSM_RNUN		dtrsm_RNUN
#define	DTRSM_RNLU		dtrsm_RNLU
#define	DTRSM_RNLN		dtrsm_RNLN
#define	DTRSM_RTUU		dtrsm_RTUU
#define	DTRSM_RTUN		dtrsm_RTUN
#define	DTRSM_RTLU		dtrsm_RTLU
#define	DTRSM_RTLN		dtrsm_RTLN
#define	DTRSM_RRUU		dtrsm_RNUU
#define	DTRSM_RRUN		dtrsm_RNUN
#define	DTRSM_RRLU		dtrsm_RNLU
#define	DTRSM_RRLN		dtrsm_RNLN
#define	DTRSM_RCUU		dtrsm_RTUU
#define	DTRSM_RCUN		dtrsm_RTUN
#define	DTRSM_RCLU		dtrsm_RTLU
#define	DTRSM_RCLN		dtrsm_RTLN

#define	DGEMM_THREAD_NN		dgemm_thread_nn
#define	DGEMM_THREAD_CN		dgemm_thread_tn
#define	DGEMM_THREAD_TN		dgemm_thread_tn
#define	DGEMM_THREAD_NC		dgemm_thread_nt
#define	DGEMM_THREAD_NT		dgemm_thread_nt
#define	DGEMM_THREAD_CC		dgemm_thread_tt
#define	DGEMM_THREAD_CT		dgemm_thread_tt
#define	DGEMM_THREAD_TC		dgemm_thread_tt
#define	DGEMM_THREAD_TT		dgemm_thread_tt
#define	DGEMM_THREAD_NR		dgemm_thread_nn
#define	DGEMM_THREAD_TR		dgemm_thread_tn
#define	DGEMM_THREAD_CR		dgemm_thread_tn
#define	DGEMM_THREAD_RN		dgemm_thread_nn
#define	DGEMM_THREAD_RT		dgemm_thread_nt
#define	DGEMM_THREAD_RC		dgemm_thread_nt
#define	DGEMM_THREAD_RR		dgemm_thread_nn

#define	DSYMM_THREAD_LU		dsymm_thread_LU
#define	DSYMM_THREAD_LL		dsymm_thread_LL
#define	DSYMM_THREAD_RU		dsymm_thread_RU
#define	DSYMM_THREAD_RL		dsymm_thread_RL

#define	DHEMM_THREAD_LU		dhemm_thread_LU
#define	DHEMM_THREAD_LL		dhemm_thread_LL
#define	DHEMM_THREAD_RU		dhemm_thread_RU
#define	DHEMM_THREAD_RL		dhemm_thread_RL

#define	DSYRK_THREAD_UN		dsyrk_thread_UN
#define	DSYRK_THREAD_UT		dsyrk_thread_UT
#define	DSYRK_THREAD_LN		dsyrk_thread_LN
#define	DSYRK_THREAD_LT		dsyrk_thread_LT
#define	DSYRK_THREAD_UR		dsyrk_thread_UN
#define	DSYRK_THREAD_UC		dsyrk_thread_UT
#define	DSYRK_THREAD_LR		dsyrk_thread_LN
#define	DSYRK_THREAD_LC		dsyrk_thread_LT

#define	DHERK_THREAD_UN		dsyrk_thread_UN
#define	DHERK_THREAD_UT		dsyrk_thread_UT
#define	DHERK_THREAD_LN		dsyrk_thread_LN
#define	DHERK_THREAD_LT		dsyrk_thread_LT
#define	DHERK_THREAD_UR		dsyrk_thread_UN
#define	DHERK_THREAD_UC		dsyrk_thread_UT
#define	DHERK_THREAD_LR		dsyrk_thread_LN
#define	DHERK_THREAD_LC		dsyrk_thread_LT

#endif
