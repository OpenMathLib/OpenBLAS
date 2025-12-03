/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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

#define xstr(s) str(s)
#define str(s) #s

int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b)
{
    IFLOAT *aoffset = a;
    IFLOAT *boffset = b;
    IFLOAT *boffset2  = b + m  * (n & ~(BLOCKSIZE - 1)) * 2;

    //fprintf(stderr, "%s m=%ld n=%ld lda=%ld\n", xstr(CNAME), m, n, lda);

    for (BLASLONG j = m; j > 0; j--) {
        IFLOAT *aoffset1 = aoffset;
        IFLOAT *boffset1 = boffset;

        aoffset += lda * 2;
        boffset += BLOCKSIZE * 2;

        for (BLASLONG i = n / BLOCKSIZE; i > 0; i--) {
            size_t vl = BLOCKSIZE;

            FLOAT_V_T v = VLEV_FLOAT(aoffset1, vl);
            VSEV_FLOAT(boffset1, v, vl);
            v = VLEV_FLOAT(aoffset1 + vl, vl);
            VSEV_FLOAT(boffset1 + vl, v, vl);

            aoffset1 += BLOCKSIZE * 2;
            boffset1 += BLOCKSIZE * m * 2;
        }

        if (n & (BLOCKSIZE - 1)) {
            size_t vl = n & (BLOCKSIZE - 1);

            FLOAT_V_T v = VLEV_FLOAT(aoffset1, vl);
            VSEV_FLOAT(boffset2, v, vl);
            v = VLEV_FLOAT(aoffset1 + vl, vl);
            VSEV_FLOAT(boffset2 + vl, v, vl);

            boffset2 += vl * 2;
        }
    }

    return 0;
}
