#ifndef COMMON_SB_H
#define COMMON_SB_H

#ifndef DYNAMIC_ARCH

#define SBDOT_K             sbdot_k
#define SBSTOBF16_K         sbstobf16_k
#define SBDTOBF16_K         sbdtobf16_k
#define SBF16TOS_K          sbf16tos_k
#define DBF16TOD_K          dbf16tod_k
#define SBGEMV_N_K          sbgemv_n
#define SBGEMV_T_K          sbgemv_t

#define	SBGEMM_ONCOPY		sbgemm_oncopy
#define	SBGEMM_OTCOPY		sbgemm_otcopy

#if SBGEMM_DEFAULT_UNROLL_M == SBGEMM_DEFAULT_UNROLL_N
#define	SBGEMM_INCOPY		sbgemm_oncopy
#define	SBGEMM_ITCOPY		sbgemm_otcopy
#else
#define	SBGEMM_INCOPY		sbgemm_incopy
#define	SBGEMM_ITCOPY		sbgemm_itcopy
#endif
#define	SBGEMM_BETA		sbgemm_beta
#define SBGEMM_KERNEL            sbgemm_kernel

#define SBGEMM_SMALL_MATRIX_PERMIT	sbgemm_small_matrix_permit
#else

#define SBDOT_K             OPENBLAS_DISPATCH(sbdot) -> sbdot_k
#define SBSTOBF16_K         OPENBLAS_DISPATCH(sbstobf16) -> sbstobf16_k
#define SBDTOBF16_K         OPENBLAS_DISPATCH(sbdtobf16) -> sbdtobf16_k
#define SBF16TOS_K          OPENBLAS_DISPATCH(sbf16tos) -> sbf16tos_k
#define DBF16TOD_K          OPENBLAS_DISPATCH(dbf16tod) -> dbf16tod_k
#define SBGEMV_N_K          OPENBLAS_DISPATCH(sbgemv) -> sbgemv_n
#define SBGEMV_T_K          OPENBLAS_DISPATCH(sbgemv) -> sbgemv_t

#define	SBGEMM_ONCOPY		OPENBLAS_DISPATCH(sbgemm) -> sbgemm_oncopy
#define	SBGEMM_OTCOPY		OPENBLAS_DISPATCH(sbgemm) -> sbgemm_otcopy
#define	SBGEMM_INCOPY		OPENBLAS_DISPATCH(sbgemm) -> sbgemm_incopy
#define	SBGEMM_ITCOPY		OPENBLAS_DISPATCH(sbgemm) -> sbgemm_itcopy
#define	SBGEMM_BETA		OPENBLAS_DISPATCH(sbgemm) -> sbgemm_beta
#define	SBGEMM_KERNEL		OPENBLAS_DISPATCH(sbgemm) -> sbgemm_kernel

#define SBGEMM_SMALL_MATRIX_PERMIT	OPENBLAS_DISPATCH(sbgemm) -> sbgemm_small_matrix_permit
#endif

#define SBGEMM_SMALL_KERNEL_BASE		OPENBLAS_DISPATCH_BASE(sbgemm)
#define SBGEMM_SMALL_KERNEL_NN		OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_nn)
#define SBGEMM_SMALL_KERNEL_NT		OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_nt)
#define SBGEMM_SMALL_KERNEL_TN		OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_tn)
#define SBGEMM_SMALL_KERNEL_TT		OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_tt)

#define SBGEMM_SMALL_KERNEL_B0_NN	OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_b0_nn)
#define SBGEMM_SMALL_KERNEL_B0_NT	OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_b0_nt)
#define SBGEMM_SMALL_KERNEL_B0_TN	OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_b0_tn)
#define SBGEMM_SMALL_KERNEL_B0_TT	OPENBLAS_DISPATCH_OFFSET(sbgemm, sbgemm_small_kernel_b0_tt)

#define	SBGEMM_NN		sbgemm_nn
#define	SBGEMM_CN		sbgemm_tn
#define	SBGEMM_TN		sbgemm_tn
#define	SBGEMM_NC		sbgemm_nt
#define	SBGEMM_NT		sbgemm_nt
#define	SBGEMM_CC		sbgemm_tt
#define	SBGEMM_CT		sbgemm_tt
#define	SBGEMM_TC		sbgemm_tt
#define	SBGEMM_TT		sbgemm_tt
#define	SBGEMM_NR		sbgemm_nn
#define	SBGEMM_TR		sbgemm_tn
#define	SBGEMM_CR		sbgemm_tn
#define	SBGEMM_RN		sbgemm_nn
#define	SBGEMM_RT		sbgemm_nt
#define	SBGEMM_RC		sbgemm_nt
#define	SBGEMM_RR		sbgemm_nn

#define	SBGEMM_THREAD_NN		sbgemm_thread_nn
#define	SBGEMM_THREAD_CN		sbgemm_thread_tn
#define	SBGEMM_THREAD_TN		sbgemm_thread_tn
#define	SBGEMM_THREAD_NC		sbgemm_thread_nt
#define	SBGEMM_THREAD_NT		sbgemm_thread_nt
#define	SBGEMM_THREAD_CC		sbgemm_thread_tt
#define	SBGEMM_THREAD_CT		sbgemm_thread_tt
#define	SBGEMM_THREAD_TC		sbgemm_thread_tt
#define	SBGEMM_THREAD_TT		sbgemm_thread_tt
#define	SBGEMM_THREAD_NR		sbgemm_thread_nn
#define	SBGEMM_THREAD_TR		sbgemm_thread_tn
#define	SBGEMM_THREAD_CR		sbgemm_thread_tn
#define	SBGEMM_THREAD_RN		sbgemm_thread_nn
#define	SBGEMM_THREAD_RT		sbgemm_thread_nt
#define	SBGEMM_THREAD_RC		sbgemm_thread_nt
#define	SBGEMM_THREAD_RR		sbgemm_thread_nn

#endif

