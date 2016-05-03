/*******************************************************************************
Copyright (c) 2016, The OpenBLAS Project
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
*******************************************************************************/

#include "common.h"
#include "macros_msa.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT * __restrict src, BLASLONG lda,
          FLOAT *  __restrict dst)
{
    BLASLONG i, j;
    FLOAT *psrc0, *psrc1, *psrc2, *psrc3, *psrc4;
    FLOAT *psrc5, *psrc6, *psrc7, *psrc8;
    FLOAT *pdst;
    v2f64 src0, src1, src2, src3, src4, src5, src6, src7;
    v2f64 src8, src9, src10, src11, src12, src13, src14, src15;
    v2f64 dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7;

    psrc0 = src;
    pdst = dst;

    for (j = (n >> 3); j--;)
    {
        psrc1 = psrc0;
        psrc2 = psrc1 + lda;
        psrc3 = psrc2 + lda;
        psrc4 = psrc3 + lda;
        psrc5 = psrc4 + lda;
        psrc6 = psrc5 + lda;
        psrc7 = psrc6 + lda;
        psrc8 = psrc7 + lda;
        psrc0 += 8 * lda;

        for (i = (m >> 3); i--;)
        {
            LD_DP2(psrc1, 2, src0, src1);
            LD_DP2(psrc2, 2, src2, src3);
            LD_DP2(psrc3, 2, src4, src5);
            LD_DP2(psrc4, 2, src6, src7);
            LD_DP2(psrc5, 2, src8, src9);
            LD_DP2(psrc6, 2, src10, src11);
            LD_DP2(psrc7, 2, src12, src13);
            LD_DP2(psrc8, 2, src14, src15);

            dst0 = (v2f64) __msa_ilvr_d((v2i64) src2, (v2i64) src0);
            dst1 = (v2f64) __msa_ilvr_d((v2i64) src6, (v2i64) src4);
            dst2 = (v2f64) __msa_ilvr_d((v2i64) src10, (v2i64) src8);
            dst3 = (v2f64) __msa_ilvr_d((v2i64) src14, (v2i64) src12);
            dst4 = (v2f64) __msa_ilvl_d((v2i64) src2, (v2i64) src0);
            dst5 = (v2f64) __msa_ilvl_d((v2i64) src6, (v2i64) src4);
            dst6 = (v2f64) __msa_ilvl_d((v2i64) src10, (v2i64) src8);
            dst7 = (v2f64) __msa_ilvl_d((v2i64) src14, (v2i64) src12);

            ST_DP8(dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7, pdst, 2);

            dst0 = (v2f64) __msa_ilvr_d((v2i64) src3, (v2i64) src1);
            dst1 = (v2f64) __msa_ilvr_d((v2i64) src7, (v2i64) src5);
            dst2 = (v2f64) __msa_ilvr_d((v2i64) src11, (v2i64) src9);
            dst3 = (v2f64) __msa_ilvr_d((v2i64) src15, (v2i64) src13);
            dst4 = (v2f64) __msa_ilvl_d((v2i64) src3, (v2i64) src1);
            dst5 = (v2f64) __msa_ilvl_d((v2i64) src7, (v2i64) src5);
            dst6 = (v2f64) __msa_ilvl_d((v2i64) src11, (v2i64) src9);
            dst7 = (v2f64) __msa_ilvl_d((v2i64) src15, (v2i64) src13);

            ST_DP8(dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7, pdst + 16,
                   2);

            LD_DP2(psrc1 + 4, 2, src0, src1);
            LD_DP2(psrc2 + 4, 2, src2, src3);
            LD_DP2(psrc3 + 4, 2, src4, src5);
            LD_DP2(psrc4 + 4, 2, src6, src7);
            LD_DP2(psrc5 + 4, 2, src8, src9);
            LD_DP2(psrc6 + 4, 2, src10, src11);
            LD_DP2(psrc7 + 4, 2, src12, src13);
            LD_DP2(psrc8 + 4, 2, src14, src15);

            dst0 = (v2f64) __msa_ilvr_d((v2i64) src2, (v2i64) src0);
            dst1 = (v2f64) __msa_ilvr_d((v2i64) src6, (v2i64) src4);
            dst2 = (v2f64) __msa_ilvr_d((v2i64) src10, (v2i64) src8);
            dst3 = (v2f64) __msa_ilvr_d((v2i64) src14, (v2i64) src12);
            dst4 = (v2f64) __msa_ilvl_d((v2i64) src2, (v2i64) src0);
            dst5 = (v2f64) __msa_ilvl_d((v2i64) src6, (v2i64) src4);
            dst6 = (v2f64) __msa_ilvl_d((v2i64) src10, (v2i64) src8);
            dst7 = (v2f64) __msa_ilvl_d((v2i64) src14, (v2i64) src12);

            ST_DP8(dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7, pdst + 32,
                   2);

            dst0 = (v2f64) __msa_ilvr_d((v2i64) src3, (v2i64) src1);
            dst1 = (v2f64) __msa_ilvr_d((v2i64) src7, (v2i64) src5);
            dst2 = (v2f64) __msa_ilvr_d((v2i64) src11, (v2i64) src9);
            dst3 = (v2f64) __msa_ilvr_d((v2i64) src15, (v2i64) src13);
            dst4 = (v2f64) __msa_ilvl_d((v2i64) src3, (v2i64) src1);
            dst5 = (v2f64) __msa_ilvl_d((v2i64) src7, (v2i64) src5);
            dst6 = (v2f64) __msa_ilvl_d((v2i64) src11, (v2i64) src9);
            dst7 = (v2f64) __msa_ilvl_d((v2i64) src15, (v2i64) src13);

            ST_DP8(dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7, pdst + 48,
                   2);

            psrc1 += 8;
            psrc2 += 8;
            psrc3 += 8;
            psrc4 += 8;
            psrc5 += 8;
            psrc6 += 8;
            psrc7 += 8;
            psrc8 += 8;
            pdst += 64;
        }

        for (i = (m & 7); i--;)
        {
            *pdst++ = *psrc1++;
            *pdst++ = *psrc2++;
            *pdst++ = *psrc3++;
            *pdst++ = *psrc4++;
            *pdst++ = *psrc5++;
            *pdst++ = *psrc6++;
            *pdst++ = *psrc7++;
            *pdst++ = *psrc8++;
        }
    }

    if (n & 4)
    {
        psrc1 = psrc0;
        psrc2 = psrc1 + lda;
        psrc3 = psrc2 + lda;
        psrc4 = psrc3 + lda;
        psrc0 += 4 * lda;

        for (i = (m >> 2); i--;)
        {
            LD_DP2(psrc1, 2, src0, src1);
            LD_DP2(psrc2, 2, src2, src3);
            LD_DP2(psrc3, 2, src4, src5);
            LD_DP2(psrc4, 2, src6, src7);
            psrc1 += 4;
            psrc2 += 4;
            psrc3 += 4;
            psrc4 += 4;

            dst0 = (v2f64) __msa_ilvr_d((v2i64) src2, (v2i64) src0);
            dst1 = (v2f64) __msa_ilvr_d((v2i64) src6, (v2i64) src4);
            dst2 = (v2f64) __msa_ilvr_d((v2i64) src3, (v2i64) src1);
            dst3 = (v2f64) __msa_ilvr_d((v2i64) src7, (v2i64) src5);

            dst4 = (v2f64) __msa_ilvl_d((v2i64) src2, (v2i64) src0);
            dst5 = (v2f64) __msa_ilvl_d((v2i64) src6, (v2i64) src4);
            dst6 = (v2f64) __msa_ilvl_d((v2i64) src3, (v2i64) src1);
            dst7 = (v2f64) __msa_ilvl_d((v2i64) src7, (v2i64) src5);

            ST_DP8(dst0, dst1, dst4, dst5, dst2, dst3, dst6, dst7, pdst, 2);
            pdst += 16;
        }

        for (i = (m & 3); i--;)
        {
            *pdst++ = *psrc1++;
            *pdst++ = *psrc2++;
            *pdst++ = *psrc3++;
            *pdst++ = *psrc4++;
        }
    }

    if (n & 2)
    {
        psrc1 = psrc0;
        psrc2 = psrc1 + lda;
        psrc0 += 2 * lda;

        for (i = (m >> 1); i--;)
        {
            src0 = LD_DP(psrc1);
            src1 = LD_DP(psrc2);
            psrc1 += 2;
            psrc2 += 2;

            dst0 = (v2f64) __msa_ilvr_d((v2i64) src1, (v2i64) src0);
            dst1 = (v2f64) __msa_ilvl_d((v2i64) src1, (v2i64) src0);

            ST_DP2(dst0, dst1, pdst, 2);
            pdst += 4;
        }

        if (m & 1)
        {
            *pdst++ = *psrc1++;
            *pdst++ = *psrc2++;
        }
    }

    if (n & 1)
    {
        psrc1 = psrc0;

        for (i = m; i--;)
        {
            *pdst++ = *psrc1++;
        }
    }

    return 0;
}
