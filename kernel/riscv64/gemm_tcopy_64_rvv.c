/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
All rights reserved.
*****************************************************************************/

#include "common.h"

/* Bitsliced tcopy for MR=64 — panels 64,32,16,8,4,2,1. */

#if !defined(DOUBLE)
#define VSETVL(n)               __riscv_vsetvl_e32m2(n)
#define FLOAT_V_T               vfloat32m2_t
#define VLEV_FLOAT              __riscv_vle32_v_f32m2
#define VSEV_FLOAT              __riscv_vse32_v_f32m2
#else
#define VSETVL(n)               __riscv_vsetvl_e64m2(n)
#define FLOAT_V_T               vfloat64m2_t
#define VLEV_FLOAT              __riscv_vle64_v_f64m2
#define VSEV_FLOAT              __riscv_vse64_v_f64m2
#endif

static void pack_panel(BLASLONG m, BLASLONG w, IFLOAT *aoffset, BLASLONG lda, IFLOAT **boffset_p)
{
    IFLOAT *boffset = *boffset_p;
    BLASLONG k;
    size_t vl;

    for (k = 0; k < m; k++) {
        BLASLONG remain = w;
        IFLOAT *src = aoffset + k * lda;
        while (remain > 0) {
            vl = VSETVL((size_t)remain);
            FLOAT_V_T v0 = VLEV_FLOAT(src, vl);
            VSEV_FLOAT(boffset, v0, vl);
            src += vl;
            boffset += vl;
            remain -= (BLASLONG)vl;
        }
    }
    *boffset_p = boffset;
}

int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b)
{
    IFLOAT *aoffset = a;
    IFLOAT *boffset = b;
    BLASLONG widths[] = {64, 32, 16, 8, 4, 2, 1};
    int wi;

    for (wi = 0; wi < 7; wi++) {
        BLASLONG w = widths[wi];
        while (n >= w) {
            pack_panel(m, w, aoffset, lda, &boffset);
            aoffset += w;
            n -= w;
        }
    }

    return 0;
}
