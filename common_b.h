/***************************************************************************
 * Copyright (c) 2025, The OpenBLAS Project
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. Neither the name of the OpenBLAS project nor the names of
 * its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * *****************************************************************************/

#ifndef COMMON_B_H
#define COMMON_B_H

#ifndef DYNAMIC_ARCH
#define BGEMV_N_K          bgemv_n
#define BGEMV_T_K          bgemv_t

#define	BSCAL_K			bscal_k

#define	BGEMM_ONCOPY		bgemm_oncopy
#define	BGEMM_OTCOPY		bgemm_otcopy

#if BGEMM_DEFAULT_UNROLL_M == BGEMM_DEFAULT_UNROLL_N
#define	BGEMM_INCOPY		bgemm_oncopy
#define	BGEMM_ITCOPY		bgemm_otcopy
#else
#define	BGEMM_INCOPY		bgemm_incopy
#define	BGEMM_ITCOPY		bgemm_itcopy
#endif

#define BGEMM_BETA bgemm_beta
#define BGEMM_KERNEL bgemm_kernel

#else
#define BGEMV_N_K          gotoblas->bgemv_n
#define BGEMV_T_K          gotoblas->bgemv_t

#define	BSCAL_K            gotoblas->bscal_k

#define BGEMM_ONCOPY gotoblas->bgemm_oncopy
#define BGEMM_OTCOPY gotoblas->bgemm_otcopy
#define BGEMM_INCOPY gotoblas->bgemm_incopy
#define BGEMM_ITCOPY gotoblas->bgemm_itcopy
#define BGEMM_BETA gotoblas->bgemm_beta
#define BGEMM_KERNEL gotoblas->bgemm_kernel

#endif

#define BGEMM_NN bgemm_nn
#define BGEMM_CN bgemm_tn
#define BGEMM_TN bgemm_tn
#define BGEMM_NC bgemm_nt
#define BGEMM_NT bgemm_nt
#define BGEMM_CC bgemm_tt
#define BGEMM_CT bgemm_tt
#define BGEMM_TC bgemm_tt
#define BGEMM_TT bgemm_tt
#define BGEMM_NR bgemm_nn
#define BGEMM_TR bgemm_tn
#define BGEMM_CR bgemm_tn
#define BGEMM_RN bgemm_nn
#define BGEMM_RT bgemm_nt
#define BGEMM_RC bgemm_nt
#define BGEMM_RR bgemm_nn

#define BGEMM_THREAD_NN bgemm_thread_nn
#define BGEMM_THREAD_CN bgemm_thread_tn
#define BGEMM_THREAD_TN bgemm_thread_tn
#define BGEMM_THREAD_NC bgemm_thread_nt
#define BGEMM_THREAD_NT bgemm_thread_nt
#define BGEMM_THREAD_CC bgemm_thread_tt
#define BGEMM_THREAD_CT bgemm_thread_tt
#define BGEMM_THREAD_TC bgemm_thread_tt
#define BGEMM_THREAD_TT bgemm_thread_tt
#define BGEMM_THREAD_NR bgemm_thread_nn
#define BGEMM_THREAD_TR bgemm_thread_tn
#define BGEMM_THREAD_CR bgemm_thread_tn
#define BGEMM_THREAD_RN bgemm_thread_nn
#define BGEMM_THREAD_RT bgemm_thread_nt
#define BGEMM_THREAD_RC bgemm_thread_nt
#define BGEMM_THREAD_RR bgemm_thread_nn
#endif
