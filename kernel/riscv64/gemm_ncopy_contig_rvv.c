/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
All rights reserved.
*****************************************************************************/

#include "common.h"

#if !defined(DOUBLE)
#define VSETVL(n)               __riscv_vsetvl_e32m1(n)
#define FLOAT_V_T               vfloat32m1_t
#define VLEV_FLOAT              __riscv_vle32_v_f32m1
#define VSEV_FLOAT              __riscv_vse32_v_f32m1
#else
#define VSETVL(n)               __riscv_vsetvl_e64m1(n)
#define FLOAT_V_T               vfloat64m1_t
#define VLEV_FLOAT              __riscv_vle64_v_f64m1
#define VSEV_FLOAT              __riscv_vse64_v_f64m1
#endif

/* Contiguous column-panel ncopy for generate_kernel a_pack=contiguous.
 * Args: m = K, n = M_rows. Single [n × m] panel (for each k: n values). */
int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b)
{
    FLOAT *b_offset = b;
    BLASLONG k;

    for (k = 0; k < m; k++) {
        BLASLONG r = 0;
        while (r < n) {
            /* scalar along n for simplicity/correctness; K is outer */
            *b_offset++ = a[r * lda + k];
            r++;
        }
    }

    return 0;
}
