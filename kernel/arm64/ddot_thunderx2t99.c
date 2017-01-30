/***************************************************************************
Copyright (c) 2017, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/


#include "common.h"

#include <arm_neon.h>

#define	N	"x0"	/* vector length */
#define	X	"x1"	/* "X" vector address */
#define	INC_X	"x2"	/* "X" stride */
#define	Y	"x3"	/* "Y" vector address */
#define	INC_Y	"x4"	/* "Y" stride */
#define J	"x5"	/* loop variable */

#define REG0	"xzr"
#define DOTF	"d0"
#define TMPX	"d16"
#define LD1VX	"{v16.d}[0]"
#define TMPY	"d24"
#define LD1VY	"{v24.d}[0]"
#define SZ	"8"

#if defined(SMP)
extern int blas_level1_thread_with_return_value(int mode, BLASLONG m, BLASLONG n,
	BLASLONG k, void *alpha, void *a, BLASLONG lda, void *b, BLASLONG ldb,
	void *c, BLASLONG ldc, int (*function)(), int nthreads);
#endif


static FLOAT ddot_compute(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	FLOAT  dot = 0.0 ;

	if ( n < 0 )  return(dot);

	__asm__ __volatile__ (
	"	mov	"N", %[N_]			\n"
	"	mov	"X", %[X_]			\n"
	"	mov	"INC_X", %[INCX_]		\n"
	"	mov	"Y", %[Y_]			\n"
	"	mov	"INC_Y", %[INCY_]		\n"
	"	fmov	"DOTF", "REG0"			\n"
	"	fmov	d1, "REG0"			\n"
	"	fmov	d2, "REG0"			\n"
	"	fmov	d3, "REG0"			\n"
	"	fmov	d4, "REG0"			\n"
	"	fmov	d5, "REG0"			\n"
	"	fmov	d6, "REG0"			\n"
	"	fmov	d7, "REG0"			\n"

	"	cmp	"N", xzr			\n"
	"	ble	.Ldot_kernel_L999		\n"

	"	cmp	"INC_X", #1			\n"
	"	bne	.Ldot_kernel_S_BEGIN		\n"
	"	cmp	"INC_Y", #1			\n"
	"	bne	.Ldot_kernel_S_BEGIN		\n"

	".Ldot_kernel_F_BEGIN:				\n"
	"	asr	"J", "N", #5			\n"
	"	cmp	"J", xzr			\n"
	"	beq	.Ldot_kernel_F1			\n"

	"	.align 5				\n"
	".Ldot_kernel_F32:				\n"
	"	ldp	q16, q17, ["X"]			\n"
	"	ldp	q24, q25, ["Y"]			\n"
	"	ldp	q18, q19, ["X", #32]		\n"
	"	ldp	q26, q27, ["Y", #32]		\n"
	"	fmla	v0.2d, v16.2d, v24.2d		\n"
	"	fmla	v1.2d, v17.2d, v25.2d		\n"
	"	ldp	q20, q21, ["X", #64]		\n"
	"	ldp	q28, q29, ["Y", #64]		\n"
	"	fmla	v2.2d, v18.2d, v26.2d		\n"
	"	fmla	v3.2d, v19.2d, v27.2d		\n"
	"	ldp	q22, q23, ["X", #96]		\n"
	"	ldp	q30, q31, ["Y", #96]		\n"
	"	add	"Y", "Y", #128			\n"
	"	add	"X", "X", #128			\n"
	"	fmla	v4.2d, v20.2d, v28.2d		\n"
	"	fmla	v5.2d, v21.2d, v29.2d		\n"
	"	PRFM	PLDL1KEEP, ["X", #896]		\n"
	"	PRFM	PLDL1KEEP, ["Y", #896]		\n"
	"	PRFM	PLDL1KEEP, ["X", #896+64]	\n"
	"	PRFM	PLDL1KEEP, ["Y", #896+64]	\n"
	"	fmla	v6.2d, v22.2d, v30.2d		\n"
	"	fmla	v7.2d, v23.2d, v31.2d		\n"

	"	ldp	q16, q17, ["X"]			\n"
	"	ldp	q24, q25, ["Y"]			\n"
	"	ldp	q18, q19, ["X", #32]		\n"
	"	ldp	q26, q27, ["Y", #32]		\n"
	"	fmla	v0.2d, v16.2d, v24.2d		\n"
	"	fmla	v1.2d, v17.2d, v25.2d		\n"
	"	ldp	q20, q21, ["X", #64]		\n"
	"	ldp	q28, q29, ["Y", #64]		\n"
	"	fmla	v2.2d, v18.2d, v26.2d		\n"
	"	fmla	v3.2d, v19.2d, v27.2d		\n"
	"	ldp	q22, q23, ["X", #96]		\n"
	"	ldp	q30, q31, ["Y", #96]		\n"
	"	add	"Y", "Y", #128			\n"
	"	add	"X", "X", #128			\n"
	"	fmla	v4.2d, v20.2d, v28.2d		\n"
	"	fmla	v5.2d, v21.2d, v29.2d		\n"
	"	PRFM	PLDL1KEEP, ["X", #896]		\n"
	"	PRFM	PLDL1KEEP, ["Y", #896]		\n"
	"	PRFM	PLDL1KEEP, ["X", #896+64]	\n"
	"	PRFM	PLDL1KEEP, ["Y", #896+64]	\n"
	"	fmla	v6.2d, v22.2d, v30.2d		\n"
	"	fmla	v7.2d, v23.2d, v31.2d		\n"

	"	subs	"J", "J", #1			\n"
	"	bne	.Ldot_kernel_F32		\n"

	"	fadd	v0.2d, v0.2d, v1.2d		\n"
	"	fadd	v2.2d, v2.2d, v3.2d		\n"
	"	fadd	v4.2d, v4.2d, v5.2d		\n"
	"	fadd	v6.2d, v6.2d, v7.2d		\n"
	"	fadd	v0.2d, v0.2d, v2.2d		\n"
	"	fadd	v4.2d, v4.2d, v6.2d		\n"
	"	fadd	v0.2d, v0.2d, v4.2d		\n"
	"	faddp	"DOTF", v0.2d			\n"

	".Ldot_kernel_F1:				\n"
	"	ands	"J", "N", #31			\n"
	"	ble	.Ldot_kernel_L999		\n"

	".Ldot_kernel_F10:				\n"
	"	ldr	"TMPX", ["X"]			\n"
	"	ldr	"TMPY", ["Y"]			\n"
	"	add	"X", "X", #"SZ"			\n"
	"	add	"Y", "Y", #"SZ"			\n"
	"	fmadd	"DOTF", "TMPX", "TMPY", "DOTF"	\n"

	"	subs	"J", "J", #1			\n"
	"	bne	.Ldot_kernel_F10		\n"

	"	b	.Ldot_kernel_L999		\n"

	".Ldot_kernel_S_BEGIN:				\n"
	"	lsl	"INC_X", "INC_X", #3		\n"
	"	lsl	"INC_Y", "INC_Y", #3		\n"
	"	asr	"J", "N", #2			\n"
	"	cmp	"J", xzr			\n"
	"	ble	.Ldot_kernel_S1			\n"

	".Ldot_kernel_S4:				\n"
	"	ld1	"LD1VX", ["X"], "INC_X"		\n"
	"	ld1	"LD1VY", ["Y"], "INC_Y"		\n"
	"	fmadd	"DOTF", "TMPX", "TMPY", "DOTF"	\n"
	"	ld1	"LD1VX", ["X"], "INC_X"		\n"
	"	ld1	"LD1VY", ["Y"], "INC_Y"		\n"
	"	fmadd	"DOTF", "TMPX", "TMPY", "DOTF"	\n"
	"	ld1	"LD1VX", ["X"], "INC_X"		\n"
	"	ld1	"LD1VY", ["Y"], "INC_Y"		\n"
	"	fmadd	"DOTF", "TMPX", "TMPY", "DOTF"	\n"
	"	ld1	"LD1VX", ["X"], "INC_X"		\n"
	"	ld1	"LD1VY", ["Y"], "INC_Y"		\n"
	"	fmadd	"DOTF", "TMPX", "TMPY", "DOTF"	\n"
	"	subs	"J", "J", #1			\n"
	"	bne	.Ldot_kernel_S4			\n"

	".Ldot_kernel_S1:				\n"
	"	ands	"J", "N", #3			\n"
	"	ble	.Ldot_kernel_L999		\n"

	".Ldot_kernel_S10:				\n"
	"	ld1	"LD1VX", ["X"], "INC_X"		\n"
	"	ld1	"LD1VY", ["Y"], "INC_Y"		\n"
	"	fmadd	"DOTF", "TMPX", "TMPY", "DOTF"	\n"
	"	subs	"J", "J", #1			\n"
	"	bne	.Ldot_kernel_S10		\n"

	".Ldot_kernel_L999:				\n"
	"	fmov	%[DOT_], "DOTF"			\n"

	: [DOT_]  "=r" (dot)		//%0
	: [N_]    "r"  (n),		//%1
	  [X_]    "r"  (x),		//%2
	  [INCX_] "r"  (inc_x),		//%3
	  [Y_]    "r"  (y),		//%4
	  [INCY_] "r"  (inc_y)		//%5
	: "cc",
	  "memory",
	  "x0", "x1", "x2", "x3", "x4", "x5",
	  "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7"
	);

	return(dot);
}

#if defined(SMP)
static int ddot_thread_function(BLASLONG n, BLASLONG dummy0,
	BLASLONG dummy1, FLOAT dummy2, FLOAT *x, BLASLONG inc_x, FLOAT *y,
	BLASLONG inc_y, FLOAT *result, BLASLONG dummy3)
{
	*result = ddot_compute(n, x, inc_x, y, inc_y);

	return 0;
}
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
#if defined(SMP)
	int nthreads;
	FLOAT dummy_alpha;
#endif
	FLOAT dot = 0.0;

#if defined(SMP)
	nthreads = num_cpu_avail(1);

	if (inc_x == 0 || inc_y == 0)
		nthreads = 1;

	if (n <= 10000)
		nthreads = 1;

	if (nthreads == 1) {
		dot = ddot_compute(n, x, inc_x, y, inc_y);
	} else {
		int mode, i;
		char result[MAX_CPU_NUMBER * sizeof(double) * 2];
		FLOAT *ptr;

		mode = BLAS_DOUBLE  | BLAS_REAL;

		blas_level1_thread_with_return_value(mode, n, 0, 0, &dummy_alpha,
				   x, inc_x, y, inc_y, result, 0,
				   ( void *)ddot_thread_function, nthreads);

		ptr = (FLOAT *)result;
		for (i = 0; i < nthreads; i++) {
			dot = dot + (*ptr);
			ptr = (FLOAT *)(((char *)ptr) + sizeof(double) * 2);
		}
	}
#else
	dot = ddot_compute(n, x, inc_x, y, inc_y);
#endif

	return dot;
}
