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

#ifndef COMMON_SH_H
#define COMMON_SH_H

#ifndef DYNAMIC_ARCH

#define SHGEMM_ONCOPY		shgemm_oncopy
#define SHGEMM_OTCOPY		shgemm_otcopy

#if SGEMM_DEFAULT_UNROLL_M == SGEMM_DEFAULT_UNROLL_N
#define SHGEMM_INCOPY		shgemm_oncopy
#define SHGEMM_ITCOPY		shgemm_otcopy
#else
#define SHGEMM_INCOPY		shgemm_incopy
#define SHGEMM_ITCOPY		shgemm_itcopy
#endif

#define SHGEMM_BETA		    shgemm_beta
#define SHGEMM_KERNEL       shgemm_kernel

#define SHGEMV_N_K      shgemv_n
#define SHGEMV_T_K      shgemv_t


#else // #DYNAMIC_ARCH

#define SHGEMM_ONCOPY		OPENBLAS_DISPATCH(shgemm) -> shgemm_oncopy
#define SHGEMM_OTCOPY		OPENBLAS_DISPATCH(shgemm) -> shgemm_otcopy
#if SGEMM_DEFAULT_UNROLL_M == SGEMM_DEFAULT_UNROLL_N
#define SHGEMM_INCOPY		OPENBLAS_DISPATCH(shgemm) -> shgemm_oncopy
#define SHGEMM_ITCOPY		OPENBLAS_DISPATCH(shgemm) -> shgemm_otcopy
#else
#define SHGEMM_INCOPY		OPENBLAS_DISPATCH(shgemm) -> shgemm_incopy
#define SHGEMM_ITCOPY		OPENBLAS_DISPATCH(shgemm) -> shgemm_itcopy
#endif

#define SHGEMM_BETA		OPENBLAS_DISPATCH(shgemm) -> shgemm_beta
#define SHGEMM_KERNEL		OPENBLAS_DISPATCH(shgemm) -> shgemm_kernel

#define SHGEMV_N_K      OPENBLAS_DISPATCH(shgemv)->shgemv_n
#define SHGEMV_T_K      OPENBLAS_DISPATCH(shgemv)->shgemv_t

#endif // #DYNAMIC_ARCH

#define SHGEMM_NN        shgemm_nn
#define SHGEMM_CN        shgemm_tn
#define SHGEMM_TN        shgemm_tn
#define SHGEMM_NC        shgemm_nt
#define SHGEMM_NT        shgemm_nt
#define SHGEMM_CC        shgemm_tt
#define SHGEMM_CT        shgemm_tt
#define SHGEMM_TC        shgemm_tt
#define SHGEMM_TT        shgemm_tt
#define SHGEMM_NR        shgemm_nn
#define SHGEMM_TR        shgemm_tn
#define SHGEMM_CR        shgemm_tn
#define SHGEMM_RN        shgemm_nn
#define SHGEMM_RT        shgemm_nt
#define SHGEMM_RC        shgemm_nt
#define SHGEMM_RR        shgemm_nn

#define SHGEMM_THREAD_NN        shgemm_thread_nn
#define SHGEMM_THREAD_CN        shgemm_thread_tn
#define SHGEMM_THREAD_TN        shgemm_thread_tn
#define SHGEMM_THREAD_NC        shgemm_thread_nt
#define SHGEMM_THREAD_NT        shgemm_thread_nt
#define SHGEMM_THREAD_CC        shgemm_thread_tt
#define SHGEMM_THREAD_CT        shgemm_thread_tt
#define SHGEMM_THREAD_TC        shgemm_thread_tt
#define SHGEMM_THREAD_TT        shgemm_thread_tt
#define SHGEMM_THREAD_NR        shgemm_thread_nn
#define SHGEMM_THREAD_TR        shgemm_thread_tn
#define SHGEMM_THREAD_CR        shgemm_thread_tn
#define SHGEMM_THREAD_RN        shgemm_thread_nn
#define SHGEMM_THREAD_RT        shgemm_thread_nt
#define SHGEMM_THREAD_RC        shgemm_thread_nt
#define SHGEMM_THREAD_RR        shgemm_thread_nn


#endif // #COMMON_SH_H