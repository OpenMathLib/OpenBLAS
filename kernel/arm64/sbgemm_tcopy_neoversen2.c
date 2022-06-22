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

#include "common.h"


int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b) {
    IFLOAT *a_offset, *a_offset1, *a_offset2, *a_offset3, *a_offset4;
    IFLOAT *b_offset;
    a_offset = a;
    b_offset = b;

    for (BLASLONG j = 0; j < n / 2; j++) {
        a_offset1 = a_offset;
        a_offset2 = a_offset1 + lda;
        a_offset3 = a_offset2 + lda;
        a_offset4 = a_offset3 + lda;
        a_offset += 2;

        for (BLASLONG i = 0; i < m / 4; i++) {
            *(b_offset + 0) = *(a_offset1 + 0);
            *(b_offset + 1) = *(a_offset2 + 0);
            *(b_offset + 2) = *(a_offset3 + 0);
            *(b_offset + 3) = *(a_offset4 + 0);
            *(b_offset + 4) = *(a_offset1 + 1);
            *(b_offset + 5) = *(a_offset2 + 1);
            *(b_offset + 6) = *(a_offset3 + 1);
            *(b_offset + 7) = *(a_offset4 + 1);

            b_offset += 8;
            a_offset1 += 4 * lda;
            a_offset2 += 4 * lda;
            a_offset3 += 4 * lda;
            a_offset4 += 4 * lda;
        }
        
        if (m & 3) {
            BLASLONG rest = m & 3;
            if (rest == 3) {
                *(b_offset + 0) = *(a_offset1 + 0);
                *(b_offset + 1) = *(a_offset2 + 0);
                *(b_offset + 2) = *(a_offset3 + 0);
                *(b_offset + 3) = *(a_offset1 + 1);
                *(b_offset + 4) = *(a_offset2 + 1);
                *(b_offset + 5) = *(a_offset3 + 1);
                b_offset += 6;
            } else if (rest == 2) {
                *(b_offset + 0) = *(a_offset1 + 0);
                *(b_offset + 1) = *(a_offset2 + 0);
                *(b_offset + 2) = *(a_offset1 + 1);
                *(b_offset + 3) = *(a_offset2 + 1);
                b_offset += 4;
            } else if (rest == 1) {
                *(b_offset + 0) = *(a_offset1 + 0);
                *(b_offset + 1) = *(a_offset1 + 1);
                b_offset += 2;
            }
        }
    }
    if (n & 1) {
        for (BLASLONG i = 0; i < m / 4; i++) {
            *(b_offset + 0) = *(a_offset);
            *(b_offset + 1) = *(a_offset + lda);
            *(b_offset + 2) = *(a_offset + lda * 2);
            *(b_offset + 3) = *(a_offset + lda * 3);

            b_offset += 4;
            a_offset += 4 * lda;
        }
        BLASLONG rest = m & 3;
        if (rest == 3) {
            *(b_offset + 0) = *(a_offset);
            *(b_offset + 1) = *(a_offset + lda);
            *(b_offset + 2) = *(a_offset + lda * 2);
        } else if (rest == 2) {
            *(b_offset + 0) = *(a_offset);
            *(b_offset + 1) = *(a_offset + lda);
        } else if (rest == 1) {
            *(b_offset + 0) = *(a_offset);
        }
    }

    return 0;
}
