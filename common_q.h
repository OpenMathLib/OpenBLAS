#ifndef COMMON_Q_H
#define COMMON_Q_H

#ifndef DYNAMIC_ARCH

#define	QAMAX_K			qamax_k
#define	QAMIN_K			qamin_k
#define	QMAX_K			qmax_k
#define	QMIN_K			qmin_k
#define	IQAMAX_K		iqamax_k
#define	IQAMIN_K		iqamin_k
#define	IQMAX_K			iqmax_k
#define	IQMIN_K			iqmin_k
#define	QASUM_K			qasum_k
#define	QAXPYU_K		qaxpy_k
#define	QAXPYC_K		qaxpy_k
#define	QCOPY_K			qcopy_k
#define	QDOTU_K			qdot_k
#define	QDOTC_K			qdot_k
#define	QNRM2_K			qnrm2_k
#define	QSCAL_K			qscal_k
#define	QSUM_K			qsum_k
#define	QSWAP_K			qswap_k
#define	QROT_K			qrot_k
#define	QROTM_K			qrotm_k

#define	QGEMV_N			qgemv_n
#define	QGEMV_T			qgemv_t
#define	QGEMV_R			qgemv_n
#define	QGEMV_C			qgemv_t
#define	QGEMV_O			qgemv_n
#define	QGEMV_U			qgemv_t
#define	QGEMV_S			qgemv_n
#define	QGEMV_D			qgemv_t

#define	QGERU_K			qger_k
#define	QGERC_K			qger_k
#define	QGERV_K			qger_k
#define	QGERD_K			qger_k

#define QSYMV_U			qsymv_U
#define QSYMV_L			qsymv_L
#define QSYMV_THREAD_U		qsymv_thread_U
#define QSYMV_THREAD_L		qsymv_thread_L

#define	QGEMM_ONCOPY		qgemm_oncopy
#define	QGEMM_OTCOPY		qgemm_otcopy

#if QGEMM_DEFAULT_UNROLL_M == QGEMM_DEFAULT_UNROLL_N
#define	QGEMM_INCOPY		qgemm_oncopy
#define	QGEMM_ITCOPY		qgemm_otcopy
#else
#define	QGEMM_INCOPY		qgemm_incopy
#define	QGEMM_ITCOPY		qgemm_itcopy
#endif

#define	QTRMM_OUNUCOPY		qtrmm_ounucopy
#define	QTRMM_OUNNCOPY		qtrmm_ounncopy
#define	QTRMM_OUTUCOPY		qtrmm_outucopy
#define	QTRMM_OUTNCOPY		qtrmm_outncopy
#define	QTRMM_OLNUCOPY		qtrmm_olnucopy
#define	QTRMM_OLNNCOPY		qtrmm_olnncopy
#define	QTRMM_OLTUCOPY		qtrmm_oltucopy
#define	QTRMM_OLTNCOPY		qtrmm_oltncopy

#define	QTRSM_OUNUCOPY		qtrsm_ounucopy
#define	QTRSM_OUNNCOPY		qtrsm_ounncopy
#define	QTRSM_OUTUCOPY		qtrsm_outucopy
#define	QTRSM_OUTNCOPY		qtrsm_outncopy
#define	QTRSM_OLNUCOPY		qtrsm_olnucopy
#define	QTRSM_OLNNCOPY		qtrsm_olnncopy
#define	QTRSM_OLTUCOPY		qtrsm_oltucopy
#define	QTRSM_OLTNCOPY		qtrsm_oltncopy

#if QGEMM_DEFAULT_UNROLL_M == QGEMM_DEFAULT_UNROLL_N
#define	QTRMM_IUNUCOPY		qtrmm_ounucopy
#define	QTRMM_IUNNCOPY		qtrmm_ounncopy
#define	QTRMM_IUTUCOPY		qtrmm_outucopy
#define	QTRMM_IUTNCOPY		qtrmm_outncopy
#define	QTRMM_ILNUCOPY		qtrmm_olnucopy
#define	QTRMM_ILNNCOPY		qtrmm_olnncopy
#define	QTRMM_ILTUCOPY		qtrmm_oltucopy
#define	QTRMM_ILTNCOPY		qtrmm_oltncopy

#define	QTRSM_IUNUCOPY		qtrsm_ounucopy
#define	QTRSM_IUNNCOPY		qtrsm_ounncopy
#define	QTRSM_IUTUCOPY		qtrsm_outucopy
#define	QTRSM_IUTNCOPY		qtrsm_outncopy
#define	QTRSM_ILNUCOPY		qtrsm_olnucopy
#define	QTRSM_ILNNCOPY		qtrsm_olnncopy
#define	QTRSM_ILTUCOPY		qtrsm_oltucopy
#define	QTRSM_ILTNCOPY		qtrsm_oltncopy
#else
#define	QTRMM_IUNUCOPY		qtrmm_iunucopy
#define	QTRMM_IUNNCOPY		qtrmm_iunncopy
#define	QTRMM_IUTUCOPY		qtrmm_iutucopy
#define	QTRMM_IUTNCOPY		qtrmm_iutncopy
#define	QTRMM_ILNUCOPY		qtrmm_ilnucopy
#define	QTRMM_ILNNCOPY		qtrmm_ilnncopy
#define	QTRMM_ILTUCOPY		qtrmm_iltucopy
#define	QTRMM_ILTNCOPY		qtrmm_iltncopy

#define	QTRSM_IUNUCOPY		qtrsm_iunucopy
#define	QTRSM_IUNNCOPY		qtrsm_iunncopy
#define	QTRSM_IUTUCOPY		qtrsm_iutucopy
#define	QTRSM_IUTNCOPY		qtrsm_iutncopy
#define	QTRSM_ILNUCOPY		qtrsm_ilnucopy
#define	QTRSM_ILNNCOPY		qtrsm_ilnncopy
#define	QTRSM_ILTUCOPY		qtrsm_iltucopy
#define	QTRSM_ILTNCOPY		qtrsm_iltncopy
#endif

#define	QGEMM_BETA		qgemm_beta

#define	QGEMM_KERNEL		qgemm_kernel

#define	QTRMM_KERNEL_LN		qtrmm_kernel_LN
#define	QTRMM_KERNEL_LT		qtrmm_kernel_LT
#define	QTRMM_KERNEL_LR		qtrmm_kernel_LN
#define	QTRMM_KERNEL_LC		qtrmm_kernel_LT
#define	QTRMM_KERNEL_RN		qtrmm_kernel_RN
#define	QTRMM_KERNEL_RT		qtrmm_kernel_RT
#define	QTRMM_KERNEL_RR		qtrmm_kernel_RN
#define	QTRMM_KERNEL_RC		qtrmm_kernel_RT

#define	QTRSM_KERNEL_LN		qtrsm_kernel_LN
#define	QTRSM_KERNEL_LT		qtrsm_kernel_LT
#define	QTRSM_KERNEL_LR		qtrsm_kernel_LN
#define	QTRSM_KERNEL_LC		qtrsm_kernel_LT
#define	QTRSM_KERNEL_RN		qtrsm_kernel_RN
#define	QTRSM_KERNEL_RT		qtrsm_kernel_RT
#define	QTRSM_KERNEL_RR		qtrsm_kernel_RN
#define	QTRSM_KERNEL_RC		qtrsm_kernel_RT

#define	QSYMM_OUTCOPY		qsymm_outcopy
#define	QSYMM_OLTCOPY		qsymm_oltcopy
#if QGEMM_DEFAULT_UNROLL_M == QGEMM_DEFAULT_UNROLL_N
#define	QSYMM_IUTCOPY		qsymm_outcopy
#define	QSYMM_ILTCOPY		qsymm_oltcopy
#else
#define	QSYMM_IUTCOPY		qsymm_iutcopy
#define	QSYMM_ILTCOPY		qsymm_iltcopy
#endif

#define QNEG_TCOPY		qneg_tcopy
#define QLASWP_NCOPY		qlaswp_ncopy

#else

#define	QAMAX_K			OPENBLAS_DISPATCH(qamax) -> qamax_k
#define	QAMIN_K			OPENBLAS_DISPATCH(qamin) -> qamin_k
#define	QMAX_K			OPENBLAS_DISPATCH(qmax) -> qmax_k
#define	QMIN_K			OPENBLAS_DISPATCH(qmin) -> qmin_k
#define	IQAMAX_K		OPENBLAS_DISPATCH(iqamax) -> iqamax_k
#define	IQAMIN_K		OPENBLAS_DISPATCH(iqamin) -> iqamin_k
#define	IQMAX_K			OPENBLAS_DISPATCH(iqmax) -> iqmax_k
#define	IQMIN_K			OPENBLAS_DISPATCH(iqmin) -> iqmin_k
#define	QASUM_K			OPENBLAS_DISPATCH(qasum) -> qasum_k
#define	QAXPYU_K		OPENBLAS_DISPATCH(qaxpy) -> qaxpy_k
#define	QAXPYC_K		OPENBLAS_DISPATCH(qaxpy) -> qaxpy_k
#define	QCOPY_K			OPENBLAS_DISPATCH(qcopy) -> qcopy_k
#define	QDOTU_K			OPENBLAS_DISPATCH(qdot) -> qdot_k
#define	QDOTC_K			OPENBLAS_DISPATCH(qdot) -> qdot_k
#define	QNRM2_K			OPENBLAS_DISPATCH(qnrm2) -> qnrm2_k
#define	QSCAL_K			OPENBLAS_DISPATCH(qscal) -> qscal_k
#define	QSUM_K			OPENBLAS_DISPATCH(qsum) -> qsum_k
#define	QSWAP_K			OPENBLAS_DISPATCH(qswap) -> qswap_k
#define	QROT_K			OPENBLAS_DISPATCH(qrot) -> qrot_k
#define	QROTM_K			OPENBLAS_DISPATCH(qrotm) -> qrotm_k

#define	QGEMV_N			OPENBLAS_DISPATCH(qgemv) -> qgemv_n
#define	QGEMV_T			OPENBLAS_DISPATCH(qgemv) -> qgemv_t
#define	QGEMV_R			OPENBLAS_DISPATCH(qgemv) -> qgemv_n
#define	QGEMV_C			OPENBLAS_DISPATCH(qgemv) -> qgemv_t
#define	QGEMV_O			OPENBLAS_DISPATCH(qgemv) -> qgemv_n
#define	QGEMV_U			OPENBLAS_DISPATCH(qgemv) -> qgemv_t
#define	QGEMV_S			OPENBLAS_DISPATCH(qgemv) -> qgemv_n
#define	QGEMV_D			OPENBLAS_DISPATCH(qgemv) -> qgemv_t

#define	QGERU_K			OPENBLAS_DISPATCH(qger) -> qger_k
#define	QGERC_K			OPENBLAS_DISPATCH(qger) -> qger_k
#define	QGERV_K			OPENBLAS_DISPATCH(qger) -> qger_k
#define	QGERD_K			OPENBLAS_DISPATCH(qger) -> qger_k

#define QSYMV_U			OPENBLAS_DISPATCH(qsymv) -> qsymv_U
#define QSYMV_L			OPENBLAS_DISPATCH(qsymv) -> qsymv_L

#define QSYMV_THREAD_U		qsymv_thread_U
#define QSYMV_THREAD_L		qsymv_thread_L

#define	QGEMM_ONCOPY		OPENBLAS_DISPATCH(qgemm) -> qgemm_oncopy
#define	QGEMM_OTCOPY		OPENBLAS_DISPATCH(qgemm) -> qgemm_otcopy
#define	QGEMM_INCOPY		OPENBLAS_DISPATCH(qgemm) -> qgemm_incopy
#define	QGEMM_ITCOPY		OPENBLAS_DISPATCH(qgemm) -> qgemm_itcopy

#define	QTRMM_OUNUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_ounucopy
#define	QTRMM_OUTUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_outucopy
#define	QTRMM_OLNUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_olnucopy
#define	QTRMM_OLTUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_oltucopy
#define	QTRSM_OUNUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_ounucopy
#define	QTRSM_OUTUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_outucopy
#define	QTRSM_OLNUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_olnucopy
#define	QTRSM_OLTUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_oltucopy

#define	QTRMM_IUNUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_iunucopy
#define	QTRMM_IUTUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_iutucopy
#define	QTRMM_ILNUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_ilnucopy
#define	QTRMM_ILTUCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_iltucopy
#define	QTRSM_IUNUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_iunucopy
#define	QTRSM_IUTUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_iutucopy
#define	QTRSM_ILNUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_ilnucopy
#define	QTRSM_ILTUCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_iltucopy

#define	QTRMM_OUNNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_ounncopy
#define	QTRMM_OUTNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_outncopy
#define	QTRMM_OLNNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_olnncopy
#define	QTRMM_OLTNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_oltncopy
#define	QTRSM_OUNNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_ounncopy
#define	QTRSM_OUTNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_outncopy
#define	QTRSM_OLNNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_olnncopy
#define	QTRSM_OLTNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_oltncopy

#define	QTRMM_IUNNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_iunncopy
#define	QTRMM_IUTNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_iutncopy
#define	QTRMM_ILNNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_ilnncopy
#define	QTRMM_ILTNCOPY		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_iltncopy
#define	QTRSM_IUNNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_iunncopy
#define	QTRSM_IUTNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_iutncopy
#define	QTRSM_ILNNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_ilnncopy
#define	QTRSM_ILTNCOPY		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_iltncopy

#define	QGEMM_BETA		OPENBLAS_DISPATCH(qgemm) -> qgemm_beta
#define	QGEMM_KERNEL		OPENBLAS_DISPATCH(qgemm) -> qgemm_kernel

#define	QTRMM_KERNEL_LN		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_LN
#define	QTRMM_KERNEL_LT		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_LT
#define	QTRMM_KERNEL_LR		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_LN
#define	QTRMM_KERNEL_LC		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_LT
#define	QTRMM_KERNEL_RN		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_RN
#define	QTRMM_KERNEL_RT		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_RT
#define	QTRMM_KERNEL_RR		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_RN
#define	QTRMM_KERNEL_RC		OPENBLAS_DISPATCH(qtrmm) -> qtrmm_kernel_RT

#define	QTRSM_KERNEL_LN		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_LN
#define	QTRSM_KERNEL_LT		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_LT
#define	QTRSM_KERNEL_LR		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_LN
#define	QTRSM_KERNEL_LC		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_LT
#define	QTRSM_KERNEL_RN		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_RN
#define	QTRSM_KERNEL_RT		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_RT
#define	QTRSM_KERNEL_RR		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_RN
#define	QTRSM_KERNEL_RC		OPENBLAS_DISPATCH(qtrsm) -> qtrsm_kernel_RT

#define	QSYMM_IUTCOPY		OPENBLAS_DISPATCH(qsymm) -> qsymm_iutcopy
#define	QSYMM_ILTCOPY		OPENBLAS_DISPATCH(qsymm) -> qsymm_iltcopy
#define	QSYMM_OUTCOPY		OPENBLAS_DISPATCH(qsymm) -> qsymm_outcopy
#define	QSYMM_OLTCOPY		OPENBLAS_DISPATCH(qsymm) -> qsymm_oltcopy

#define QNEG_TCOPY		OPENBLAS_DISPATCH(qneg) -> qneg_tcopy
#define QLASWP_NCOPY		OPENBLAS_DISPATCH(qlaswp) -> qlaswp_ncopy

#endif

#define	QGEMM_NN		qgemm_nn
#define	QGEMM_CN		qgemm_tn
#define	QGEMM_TN		qgemm_tn
#define	QGEMM_NC		qgemm_nt
#define	QGEMM_NT		qgemm_nt
#define	QGEMM_CC		qgemm_tt
#define	QGEMM_CT		qgemm_tt
#define	QGEMM_TC		qgemm_tt
#define	QGEMM_TT		qgemm_tt
#define	QGEMM_NR		qgemm_nn
#define	QGEMM_TR		qgemm_tn
#define	QGEMM_CR		qgemm_tn
#define	QGEMM_RN		qgemm_nn
#define	QGEMM_RT		qgemm_nt
#define	QGEMM_RC		qgemm_nt
#define	QGEMM_RR		qgemm_nn

#define	QSYMM_LU		qsymm_LU
#define	QSYMM_LL		qsymm_LL
#define	QSYMM_RU		qsymm_RU
#define	QSYMM_RL		qsymm_RL

#define	QHEMM_LU		qhemm_LU
#define	QHEMM_LL		qhemm_LL
#define	QHEMM_RU		qhemm_RU
#define	QHEMM_RL		qhemm_RL

#define	QSYRK_UN		qsyrk_UN
#define	QSYRK_UT		qsyrk_UT
#define	QSYRK_LN		qsyrk_LN
#define	QSYRK_LT		qsyrk_LT
#define	QSYRK_UR		qsyrk_UN
#define	QSYRK_UC		qsyrk_UT
#define	QSYRK_LR		qsyrk_LN
#define	QSYRK_LC		qsyrk_LT

#define	QSYRK_KERNEL_U		qsyrk_kernel_U
#define	QSYRK_KERNEL_L		qsyrk_kernel_L

#define	QGEMMT_UNN		qgemmt_UNN
#define	QGEMMT_UNT		qgemmt_UNT
#define	QGEMMT_UTN		qgemmt_UTN
#define	QGEMMT_UTT		qgemmt_UTT
#define	QGEMMT_LNN		qgemmt_LNN
#define	QGEMMT_LNT		qgemmt_LNT
#define	QGEMMT_LTN		qgemmt_LTN
#define	QGEMMT_LTT		qgemmt_LTT

#define	QHERK_UN		qsyrk_UN
#define	QHERK_LN		qsyrk_LN
#define	QHERK_UC		qsyrk_UT
#define	QHERK_LC		qsyrk_LT

#define	QHER2K_UN		qsyr2k_UN
#define	QHER2K_LN		qsyr2k_LN
#define	QHER2K_UC		qsyr2k_UT
#define	QHER2K_LC		qsyr2k_LT

#define	QSYR2K_UN		qsyr2k_UN
#define	QSYR2K_UT		qsyr2k_UT
#define	QSYR2K_LN		qsyr2k_LN
#define	QSYR2K_LT		qsyr2k_LT
#define	QSYR2K_UR		qsyr2k_UN
#define	QSYR2K_UC		qsyr2k_UT
#define	QSYR2K_LR		qsyr2k_LN
#define	QSYR2K_LC		qsyr2k_LT

#define	QSYR2K_KERNEL_U		qsyr2k_kernel_U
#define	QSYR2K_KERNEL_L		qsyr2k_kernel_L

#define	QTRMM_LNUU		qtrmm_LNUU
#define	QTRMM_LNUN		qtrmm_LNUN
#define	QTRMM_LNLU		qtrmm_LNLU
#define	QTRMM_LNLN		qtrmm_LNLN
#define	QTRMM_LTUU		qtrmm_LTUU
#define	QTRMM_LTUN		qtrmm_LTUN
#define	QTRMM_LTLU		qtrmm_LTLU
#define	QTRMM_LTLN		qtrmm_LTLN
#define	QTRMM_LRUU		qtrmm_LNUU
#define	QTRMM_LRUN		qtrmm_LNUN
#define	QTRMM_LRLU		qtrmm_LNLU
#define	QTRMM_LRLN		qtrmm_LNLN
#define	QTRMM_LCUU		qtrmm_LTUU
#define	QTRMM_LCUN		qtrmm_LTUN
#define	QTRMM_LCLU		qtrmm_LTLU
#define	QTRMM_LCLN		qtrmm_LTLN
#define	QTRMM_RNUU		qtrmm_RNUU
#define	QTRMM_RNUN		qtrmm_RNUN
#define	QTRMM_RNLU		qtrmm_RNLU
#define	QTRMM_RNLN		qtrmm_RNLN
#define	QTRMM_RTUU		qtrmm_RTUU
#define	QTRMM_RTUN		qtrmm_RTUN
#define	QTRMM_RTLU		qtrmm_RTLU
#define	QTRMM_RTLN		qtrmm_RTLN
#define	QTRMM_RRUU		qtrmm_RNUU
#define	QTRMM_RRUN		qtrmm_RNUN
#define	QTRMM_RRLU		qtrmm_RNLU
#define	QTRMM_RRLN		qtrmm_RNLN
#define	QTRMM_RCUU		qtrmm_RTUU
#define	QTRMM_RCUN		qtrmm_RTUN
#define	QTRMM_RCLU		qtrmm_RTLU
#define	QTRMM_RCLN		qtrmm_RTLN

#define	QTRSM_LNUU		qtrsm_LNUU
#define	QTRSM_LNUN		qtrsm_LNUN
#define	QTRSM_LNLU		qtrsm_LNLU
#define	QTRSM_LNLN		qtrsm_LNLN
#define	QTRSM_LTUU		qtrsm_LTUU
#define	QTRSM_LTUN		qtrsm_LTUN
#define	QTRSM_LTLU		qtrsm_LTLU
#define	QTRSM_LTLN		qtrsm_LTLN
#define	QTRSM_LRUU		qtrsm_LNUU
#define	QTRSM_LRUN		qtrsm_LNUN
#define	QTRSM_LRLU		qtrsm_LNLU
#define	QTRSM_LRLN		qtrsm_LNLN
#define	QTRSM_LCUU		qtrsm_LTUU
#define	QTRSM_LCUN		qtrsm_LTUN
#define	QTRSM_LCLU		qtrsm_LTLU
#define	QTRSM_LCLN		qtrsm_LTLN
#define	QTRSM_RNUU		qtrsm_RNUU
#define	QTRSM_RNUN		qtrsm_RNUN
#define	QTRSM_RNLU		qtrsm_RNLU
#define	QTRSM_RNLN		qtrsm_RNLN
#define	QTRSM_RTUU		qtrsm_RTUU
#define	QTRSM_RTUN		qtrsm_RTUN
#define	QTRSM_RTLU		qtrsm_RTLU
#define	QTRSM_RTLN		qtrsm_RTLN
#define	QTRSM_RRUU		qtrsm_RNUU
#define	QTRSM_RRUN		qtrsm_RNUN
#define	QTRSM_RRLU		qtrsm_RNLU
#define	QTRSM_RRLN		qtrsm_RNLN
#define	QTRSM_RCUU		qtrsm_RTUU
#define	QTRSM_RCUN		qtrsm_RTUN
#define	QTRSM_RCLU		qtrsm_RTLU
#define	QTRSM_RCLN		qtrsm_RTLN

#define	QGEMM_THREAD_NN		qgemm_thread_nn
#define	QGEMM_THREAD_CN		qgemm_thread_tn
#define	QGEMM_THREAD_TN		qgemm_thread_tn
#define	QGEMM_THREAD_NC		qgemm_thread_nt
#define	QGEMM_THREAD_NT		qgemm_thread_nt
#define	QGEMM_THREAD_CC		qgemm_thread_tt
#define	QGEMM_THREAD_CT		qgemm_thread_tt
#define	QGEMM_THREAD_TC		qgemm_thread_tt
#define	QGEMM_THREAD_TT		qgemm_thread_tt
#define	QGEMM_THREAD_NR		qgemm_thread_nn
#define	QGEMM_THREAD_TR		qgemm_thread_tn
#define	QGEMM_THREAD_CR		qgemm_thread_tn
#define	QGEMM_THREAD_RN		qgemm_thread_nn
#define	QGEMM_THREAD_RT		qgemm_thread_nt
#define	QGEMM_THREAD_RC		qgemm_thread_nt
#define	QGEMM_THREAD_RR		qgemm_thread_nn

#define	QSYMM_THREAD_LU		qsymm_thread_LU
#define	QSYMM_THREAD_LL		qsymm_thread_LL
#define	QSYMM_THREAD_RU		qsymm_thread_RU
#define	QSYMM_THREAD_RL		qsymm_thread_RL

#define	QHEMM_THREAD_LU		qhemm_thread_LU
#define	QHEMM_THREAD_LL		qhemm_thread_LL
#define	QHEMM_THREAD_RU		qhemm_thread_RU
#define	QHEMM_THREAD_RL		qhemm_thread_RL

#define	QSYRK_THREAD_UN		qsyrk_thread_UN
#define	QSYRK_THREAD_UT		qsyrk_thread_UT
#define	QSYRK_THREAD_LN		qsyrk_thread_LN
#define	QSYRK_THREAD_LT		qsyrk_thread_LT
#define	QSYRK_THREAD_UR		qsyrk_thread_UN
#define	QSYRK_THREAD_UC		qsyrk_thread_UT
#define	QSYRK_THREAD_LR		qsyrk_thread_LN
#define	QSYRK_THREAD_LC		qsyrk_thread_LT

#define	QHERK_THREAD_UN		qsyrk_thread_UN
#define	QHERK_THREAD_UT		qsyrk_thread_UT
#define	QHERK_THREAD_LN		qsyrk_thread_LN
#define	QHERK_THREAD_LT		qsyrk_thread_LT
#define	QHERK_THREAD_UR		qsyrk_thread_UN
#define	QHERK_THREAD_UC		qsyrk_thread_UT
#define	QHERK_THREAD_LR		qsyrk_thread_LN
#define	QHERK_THREAD_LC		qsyrk_thread_LT

#endif
