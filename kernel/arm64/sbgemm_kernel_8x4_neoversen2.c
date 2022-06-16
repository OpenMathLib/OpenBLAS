/***************************************************************************
 * Copyright (c) 2022, The OpenBLAS Project
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

#include <arm_sve.h>

#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, IFLOAT *A, IFLOAT *B, FLOAT *C,
          BLASLONG ldc) {
    // printf("m: %d, n: %d, k: %d\n", m, n, k);
    BLASLONG padk = (k + 3) & ~3;
    BLASLONG padm = (m + 1) & ~1;
    BLASLONG padn = (n + 1) & ~1;
    FLOAT *RC = (FLOAT *)calloc(padm * padn, sizeof(float));
    BLASLONG nldc = padm;

    IFLOAT *ptr_a = A;
    IFLOAT *ptr_b = B;
    FLOAT *ptr_c = RC;

    IFLOAT *ptr_a0, *ptr_a1, *ptr_a2, *ptr_a3;
    IFLOAT *ptr_b0, *ptr_b1;
    FLOAT *ptr_c00, *ptr_c10, *ptr_c20, *ptr_c30, *ptr_c01, *ptr_c11, *ptr_c21, *ptr_c31;

    svbfloat16_t ma0, ma1, ma2, ma3, mb0, mb1;
    svfloat32_t mc00, mc01, mc10, mc11, mc20, mc21, mc30, mc31;
    svbool_t pg16 = svptrue_b16();
    svbool_t pg32 = svptrue_b32();
    svfloat32_t svalpha = svdup_f32(alpha);

    uint32_t off_c[] = {0, (uint32_t)nldc, 1, (uint32_t)nldc + 1};  // 00 01 10 11
    svuint32_t off_vc = svld1_u32(pg32, off_c);

    for (BLASLONG j = 0; j < padn / 4; j++) {
        ptr_c00 = ptr_c;
        ptr_c10 = ptr_c00 + 2;
        ptr_c20 = ptr_c10 + 2;
        ptr_c30 = ptr_c20 + 2;
        ptr_c01 = ptr_c + 2 * nldc;
        ptr_c11 = ptr_c01 + 2;
        ptr_c21 = ptr_c11 + 2;
        ptr_c31 = ptr_c21 + 2;
        ptr_c += 4 * nldc;

        ptr_a = A;

        for (BLASLONG i = 0; i < padm / 8; i++) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * padk;
            ptr_a2 = ptr_a1 + 2 * padk;
            ptr_a3 = ptr_a2 + 2 * padk;
            ptr_a += 8 * padk;

            ptr_b0 = ptr_b;
            ptr_b1 = ptr_b0 + 2 * padk;

            mc00 = svdup_f32(0);
            mc01 = svdup_f32(0);
            mc10 = svdup_f32(0);
            mc11 = svdup_f32(0);
            mc20 = svdup_f32(0);
            mc21 = svdup_f32(0);
            mc30 = svdup_f32(0);
            mc31 = svdup_f32(0);

            for (BLASLONG p = 0; p < padk / 4; p++) {
                ma0 = svld1_bf16(pg16, (bfloat16_t *)ptr_a0);
                ma1 = svld1_bf16(pg16, (bfloat16_t *)ptr_a1);
                ma2 = svld1_bf16(pg16, (bfloat16_t *)ptr_a2);
                ma3 = svld1_bf16(pg16, (bfloat16_t *)ptr_a3);
                mb0 = svld1_bf16(pg16, (bfloat16_t *)ptr_b0);
                mb1 = svld1_bf16(pg16, (bfloat16_t *)ptr_b1);

                mc00 = svbfmmla(mc00, ma0, mb0);
                mc10 = svbfmmla(mc10, ma1, mb0);
                mc20 = svbfmmla(mc20, ma2, mb0);
                mc30 = svbfmmla(mc30, ma3, mb0);
                mc01 = svbfmmla(mc01, ma0, mb1);
                mc11 = svbfmmla(mc11, ma1, mb1);
                mc21 = svbfmmla(mc21, ma2, mb1);
                mc31 = svbfmmla(mc31, ma3, mb1);

                ptr_a0 += 8;
                ptr_a1 += 8;
                ptr_a2 += 8;
                ptr_a3 += 8;
                ptr_b0 += 8;
                ptr_b1 += 8;
            }
            svst1_scatter_index(pg32, ptr_c00, off_vc, mc00);
            svst1_scatter_index(pg32, ptr_c10, off_vc, mc10);
            svst1_scatter_index(pg32, ptr_c20, off_vc, mc20);
            svst1_scatter_index(pg32, ptr_c30, off_vc, mc30);
            svst1_scatter_index(pg32, ptr_c01, off_vc, mc01);
            svst1_scatter_index(pg32, ptr_c11, off_vc, mc11);
            svst1_scatter_index(pg32, ptr_c21, off_vc, mc21);
            svst1_scatter_index(pg32, ptr_c31, off_vc, mc31);

            ptr_c00 += 8;
            ptr_c10 += 8;
            ptr_c20 += 8;
            ptr_c30 += 8;
            ptr_c01 += 8;
            ptr_c11 += 8;
            ptr_c21 += 8;
            ptr_c31 += 8;
        }

        if (padm & 4) {
            // rest 4 or 6
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * padk;
            ptr_a += 4 * padk;

            ptr_b0 = ptr_b;
            ptr_b1 = ptr_b0 + 2 * padk;

            mc00 = svdup_f32(0);
            mc01 = svdup_f32(0);
            mc10 = svdup_f32(0);
            mc11 = svdup_f32(0);
            for (BLASLONG p = 0; p < padk / 4; p++) {
                ma0 = svld1_bf16(pg16, (bfloat16_t *)ptr_a0);
                ma1 = svld1_bf16(pg16, (bfloat16_t *)ptr_a1);
                mb0 = svld1_bf16(pg16, (bfloat16_t *)ptr_b0);
                mb1 = svld1_bf16(pg16, (bfloat16_t *)ptr_b1);

                mc00 = svbfmmla(mc00, ma0, mb0);
                mc10 = svbfmmla(mc10, ma1, mb0);
                mc01 = svbfmmla(mc01, ma0, mb1);
                mc11 = svbfmmla(mc11, ma1, mb1);

                ptr_a0 += 8;
                ptr_a1 += 8;
                ptr_b0 += 8;
                ptr_b1 += 8;
            }
            svst1_scatter_index(pg32, ptr_c00, off_vc, mc00);
            svst1_scatter_index(pg32, ptr_c10, off_vc, mc10);
            svst1_scatter_index(pg32, ptr_c01, off_vc, mc01);
            svst1_scatter_index(pg32, ptr_c11, off_vc, mc11);

            ptr_c00 += 4;
            ptr_c10 += 4;
            ptr_c01 += 4;
            ptr_c11 += 4;
        }

        if (padm & 2) {
            // rest 2
            ptr_a0 = ptr_a;

            ptr_b0 = ptr_b;
            ptr_b1 = ptr_b0 + 2 * padk;

            mc00 = svdup_f32(0);
            mc01 = svdup_f32(0);
            for (BLASLONG p = 0; p < padk / 4; p++) {
                ma0 = svld1_bf16(pg16, (bfloat16_t *)ptr_a0);
                mb0 = svld1_bf16(pg16, (bfloat16_t *)ptr_b0);
                mb1 = svld1_bf16(pg16, (bfloat16_t *)ptr_b1);
                mc00 = svbfmmla(mc00, ma0, mb0);
                mc01 = svbfmmla(mc01, ma0, mb1);
                ptr_a0 += 8;
                ptr_b0 += 8;
                ptr_b1 += 8;
            }
            svst1_scatter_index(pg32, ptr_c00, off_vc, mc00);
            svst1_scatter_index(pg32, ptr_c01, off_vc, mc01);
            ptr_c00 += 2;
            ptr_c01 += 2;
        }

        ptr_b += 4 * padk;
    }

    if (padn & 2) {
        // rest 2
        ptr_c00 = ptr_c;
        ptr_c10 = ptr_c00 + 2;
        ptr_c20 = ptr_c10 + 2;
        ptr_c30 = ptr_c20 + 2;
        ptr_c += 2 * nldc;

        ptr_a = A;

        for (BLASLONG i = 0; i < padm / 8; i++) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * padk;
            ptr_a2 = ptr_a1 + 2 * padk;
            ptr_a3 = ptr_a2 + 2 * padk;
            ptr_a += 8 * padk;

            ptr_b0 = ptr_b;

            mc00 = svdup_f32(0);
            mc10 = svdup_f32(0);
            mc20 = svdup_f32(0);
            mc30 = svdup_f32(0);

            for (BLASLONG p = 0; p < padk / 4; p++) {
                ma0 = svld1_bf16(pg16, (bfloat16_t *)ptr_a0);
                ma1 = svld1_bf16(pg16, (bfloat16_t *)ptr_a1);
                ma2 = svld1_bf16(pg16, (bfloat16_t *)ptr_a2);
                ma3 = svld1_bf16(pg16, (bfloat16_t *)ptr_a3);
                mb0 = svld1_bf16(pg16, (bfloat16_t *)ptr_b0);
                mc00 = svbfmmla(mc00, ma0, mb0);
                mc10 = svbfmmla(mc10, ma1, mb0);
                mc20 = svbfmmla(mc20, ma2, mb0);
                mc30 = svbfmmla(mc30, ma3, mb0);
                ptr_a0 += 8;
                ptr_a1 += 8;
                ptr_a2 += 8;
                ptr_a3 += 8;
                ptr_b0 += 8;
            }
            svst1_scatter_index(pg32, ptr_c00, off_vc, mc00);
            svst1_scatter_index(pg32, ptr_c10, off_vc, mc10);
            svst1_scatter_index(pg32, ptr_c20, off_vc, mc20);
            svst1_scatter_index(pg32, ptr_c30, off_vc, mc30);
            ptr_c00 += 8;
            ptr_c10 += 8;
            ptr_c20 += 8;
            ptr_c30 += 8;
        }

        if (padm & 4) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * padk;
            ptr_a += 4 * padk;

            ptr_b0 = ptr_b;

            mc00 = svdup_f32(0);
            mc10 = svdup_f32(0);
            for (BLASLONG p = 0; p < padk / 4; p++) {
                ma0 = svld1_bf16(pg16, (bfloat16_t *)ptr_a0);
                ma1 = svld1_bf16(pg16, (bfloat16_t *)ptr_a1);
                mb0 = svld1_bf16(pg16, (bfloat16_t *)ptr_b0);
                mc00 = svbfmmla(mc00, ma0, mb0);
                mc10 = svbfmmla(mc10, ma1, mb0);
                ptr_a0 += 8;
                ptr_a1 += 8;
                ptr_b0 += 8;
            }
            svst1_scatter_index(pg32, ptr_c00, off_vc, mc00);
            svst1_scatter_index(pg32, ptr_c10, off_vc, mc10);
            ptr_c00 += 4;
            ptr_c10 += 4;
        }

        if (padm & 2) {
            ptr_a0 = ptr_a;
            ptr_a += 2 * padk;
            ptr_b0 = ptr_b;
            mc00 = svdup_f32(0);
            for (BLASLONG p = 0; p < padk / 4; p++) {
                ma0 = svld1_bf16(pg16, (bfloat16_t *)ptr_a0);
                mb0 = svld1_bf16(pg16, (bfloat16_t *)ptr_b0);
                mc00 = svbfmmla(mc00, ma0, mb0);
                ptr_a0 += 8;
                ptr_b0 += 8;
            }
            svst1_scatter_index(pg32, ptr_c00, off_vc, mc00);
            ptr_c00 += 2;
        }

        ptr_b += 2 * padk;
    }

    FLOAT *org_c = C;
    FLOAT *raw_c = RC;
    FLOAT *org_c0, *raw_c0;
    svfloat32_t org_vc0, raw_vc0;
    for (BLASLONG j = 0; j < n; j++) {
        org_c0 = org_c;
        raw_c0 = raw_c;
        org_c += ldc;
        raw_c += nldc;
        BLASLONG i;
        for (i = 0; i < m / 4; i++) {
            org_vc0 = svld1_f32(pg32, org_c0);
            raw_vc0 = svld1_f32(pg32, raw_c0);
            org_vc0 = svmad_z(pg32, svalpha, raw_vc0,
                              org_vc0);  // alpha * raw + org, raw -> a * b
            svst1_f32(pg32, org_c0, org_vc0);
            org_c0 += 4;
            raw_c0 += 4;
        }
        for (i = 0; i < (m & 3); i++) {
            *org_c0 += alpha * (*raw_c0);
            org_c0++;
            raw_c0++;
        }
    }

    return 0;
}
