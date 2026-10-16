/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
All rights reserved.
*****************************************************************************/

#include "common.h"

/* Bitsliced symm_ucopy: pack n in panels 64,16,8,4,2,1 (a_pack=bitsliced). */

#if !defined(DOUBLE)
#define VSETVL(n)               __riscv_vsetvl_e32m2(n)
#define FLOAT_V_T               vfloat32m2_t
#define VLEV_FLOAT              __riscv_vle32_v_f32m2
#define VSEV_FLOAT              __riscv_vse32_v_f32m2
#define VLSEV_FLOAT             __riscv_vlse32_v_f32m2
#define INT_V_T                 vint32m2_t
#define VID_V_INT               __riscv_vid_v_u32m2
#define VADD_VX_INT             __riscv_vadd_vx_i32m2
#define VMSGT_VX_INT            __riscv_vmsgt_vx_i32m2_b16
#define VBOOL_T                 vbool16_t
#define VMERGE_VVM_FLOAT        __riscv_vmerge_vvm_f32m2
#define V_UM2_TO_IM2            __riscv_vreinterpret_v_u32m2_i32m2
#else
#define VSETVL(n)               __riscv_vsetvl_e64m2(n)
#define FLOAT_V_T               vfloat64m2_t
#define VLEV_FLOAT              __riscv_vle64_v_f64m2
#define VSEV_FLOAT              __riscv_vse64_v_f64m2
#define VLSEV_FLOAT             __riscv_vlse64_v_f64m2
#define INT_V_T                 vint64m2_t
#define VID_V_INT               __riscv_vid_v_u64m2
#define VADD_VX_INT             __riscv_vadd_vx_i64m2
#define VMSGT_VX_INT            __riscv_vmsgt_vx_i64m2_b32
#define VBOOL_T                 vbool32_t
#define VMERGE_VVM_FLOAT        __riscv_vmerge_vvm_f64m2
#define V_UM2_TO_IM2            __riscv_vreinterpret_v_u64m2_i64m2
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b)
{
    BLASLONG i, offset;
    FLOAT *ao1, *ao2;
    BLASLONG stride_lda = sizeof(FLOAT) * lda;
    FLOAT_V_T vb, va1, va2;
    VBOOL_T vbool;
    INT_V_T vindex_max, vindex;
    size_t vl;
    BLASLONG widths[] = {64, 32, 16, 8, 4, 2, 1};
    int wi;

    for (wi = 0; wi < 7; wi++) {
        BLASLONG w = widths[wi];
        while (n >= w) {
            vl = VSETVL((size_t)w);
            vindex_max = V_UM2_TO_IM2(VID_V_INT(vl));
            offset = posX - posY;
            ao1 = a + posY + (posX + 0) * lda;
            ao2 = a + posX + 0 + posY * lda;

            for (i = m; i > 0; i--, offset--) {
                va1 = VLSEV_FLOAT(ao1, stride_lda, vl);
                va2 = VLEV_FLOAT(ao2, vl);
                vindex = VADD_VX_INT(vindex_max, offset, vl);
                vbool = VMSGT_VX_INT(vindex, 0, vl);
                vb = VMERGE_VVM_FLOAT(va2, va1, vbool, vl);
                VSEV_FLOAT(b, vb, vl);
                b += vl;
                ao1++;
                ao2 += lda;
            }
            posX += w;
            n -= w;
        }
    }
    return 0;
}
