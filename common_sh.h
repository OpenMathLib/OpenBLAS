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


#else // #DYNAMIC_ARCH

#define SHGEMM_ONCOPY		gotoblas -> shgemm_oncopy
#define SHGEMM_OTCOPY		gotoblas -> shgemm_otcopy
#if SGEMM_DEFAULT_UNROLL_M == SGEMM_DEFAULT_UNROLL_N
#define SHGEMM_INCOPY		gotoblas -> shgemm_oncopy
#define SHGEMM_ITCOPY		gotoblas -> shgemm_otcopy
#else
#define SHGEMM_INCOPY		gotoblas -> shgemm_incopy
#define SHGEMM_ITCOPY		gotoblas -> shgemm_itcopy
#endif

#define SHGEMM_BETA		gotoblas -> shgemm_beta
#define SHGEMM_KERNEL		gotoblas -> shgemm_kernel
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