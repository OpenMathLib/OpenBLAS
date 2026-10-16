/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
All rights reserved.
*****************************************************************************/

#include "common.h"

/* Bitsliced ncopy for MR=32 — power-of-two panels matching a_pack=bitsliced. */

#if !defined(DOUBLE)
#define VSETVL(n)               __riscv_vsetvl_e32m2(n)
#define FLOAT_V_T               vfloat32m2_t
#define VLSEV_FLOAT             __riscv_vlse32_v_f32m2
#define VSEV_FLOAT              __riscv_vse32_v_f32m2
#else
#define VSETVL(n)               __riscv_vsetvl_e64m2(n)
#define FLOAT_V_T               vfloat64m2_t
#define VLSEV_FLOAT             __riscv_vlse64_v_f64m2
#define VSEV_FLOAT              __riscv_vse64_v_f64m2
#endif

static void pack_panel(BLASLONG m, BLASLONG w, FLOAT *a_offset, BLASLONG lda, FLOAT **b_offset_p)
{
    FLOAT *b_offset = *b_offset_p;
    BLASLONG k;
    size_t vl;

    for (k = 0; k < m; k++) {
        BLASLONG remain = w;
        FLOAT *src = a_offset + k;
        while (remain > 0) {
            vl = VSETVL((size_t)remain);
            FLOAT_V_T v0 = VLSEV_FLOAT(src, lda * sizeof(FLOAT), vl);
            VSEV_FLOAT(b_offset, v0, vl);
            src += (BLASLONG)vl * lda;
            b_offset += vl;
            remain -= (BLASLONG)vl;
        }
    }
    *b_offset_p = b_offset;
}

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b)
{
    /* m = K, n = M_rows. Gather rows via stride lda. */
    FLOAT *a_offset = a;
    FLOAT *b_offset = b;
    BLASLONG widths[] = {32, 16, 8, 4, 2, 1};
    int wi;

    for (wi = 0; wi < 6; wi++) {
        BLASLONG w = widths[wi];
        while (n >= w) {
            pack_panel(m, w, a_offset, lda, &b_offset);
            a_offset += w * lda;
            n -= w;
        }
    }

    return 0;
}
