/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
All rights reserved.
*****************************************************************************/

#include "common.h"

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

/* Contiguous column-panel tcopy for generate_kernel a_pack=contiguous.
 * Unlike gemm_tcopy_rvv_v1, does NOT split into VLMAX-sized sub-panels:
 * layout is a single [n × m] panel (for each k: n consecutive rows). */
int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b)
{
    IFLOAT *aoffset = a;
    IFLOAT *boffset = b;
    BLASLONG i;

    for (i = 0; i < m; i++) {
        IFLOAT *aoffset1 = aoffset;
        BLASLONG remain = n;
        while (remain > 0) {
            size_t vl = VSETVL(remain);
            FLOAT_V_T v0 = VLEV_FLOAT(aoffset1, vl);
            VSEV_FLOAT(boffset, v0, vl);
            aoffset1 += vl;
            boffset += vl;
            remain -= (BLASLONG)vl;
        }
        aoffset += lda;
    }

    return 0;
}
