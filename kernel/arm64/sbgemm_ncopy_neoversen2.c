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
    IFLOAT *a_offset, *a_offset1, *a_offset2;
    IFLOAT *b_offset;

    a_offset = a;
    b_offset = b;

    BLASLONG m4 = m & ~3;
    BLASLONG n2 = n & ~1;

    BLASLONG j = 0;
    for (; j < n2; j += 2) {
        a_offset1 = a_offset;
        a_offset2 = a_offset1 + lda;
        a_offset += 2 * lda;

        BLASLONG i = 0;
        for (; i < m4; i += 4) {
            *(b_offset + 0) = *(a_offset1 + 0);
            *(b_offset + 1) = *(a_offset1 + 1);
            *(b_offset + 2) = *(a_offset1 + 2);
            *(b_offset + 3) = *(a_offset1 + 3);
            *(b_offset + 4) = *(a_offset2 + 0);
            *(b_offset + 5) = *(a_offset2 + 1);
            *(b_offset + 6) = *(a_offset2 + 2);
            *(b_offset + 7) = *(a_offset2 + 3);

            a_offset1 += 4;
            a_offset2 += 4;
            b_offset += 8;
        }
        if (i < m) {
            *(b_offset + 0) = *(a_offset1 + 0);
            *(b_offset + 4) = *(a_offset2 + 0);

            if (i + 1 < m) {
                *(b_offset + 1) = *(a_offset1 + 1);
                *(b_offset + 5) = *(a_offset2 + 1);
            } else {
                *(b_offset + 1) = 0;
                *(b_offset + 5) = 0;
            }

            if (i + 2 < m) {
                *(b_offset + 2) = *(a_offset1 + 2);
                *(b_offset + 6) = *(a_offset2 + 2);
            } else {
                *(b_offset + 2) = 0;
                *(b_offset + 6) = 0;
            }

            *(b_offset + 3) = 0;
            *(b_offset + 7) = 0;

            b_offset += 8;
        }
    }
    if (j < n) {
        BLASLONG i = 0;
        for (; i < m4; i += 4) {
            *(b_offset + 0) = *(a_offset + 0);
            *(b_offset + 1) = *(a_offset + 1);
            *(b_offset + 2) = *(a_offset + 2);
            *(b_offset + 3) = *(a_offset + 3);
            *(b_offset + 4) = 0;
            *(b_offset + 5) = 0;
            *(b_offset + 6) = 0;
            *(b_offset + 7) = 0;
            a_offset += 4;
            b_offset += 8;
        }
        if (i < m) {
            *(b_offset + 4) = 0;
            *(b_offset + 5) = 0;
            *(b_offset + 6) = 0;
            *(b_offset + 7) = 0;

            *(b_offset + 0) = *(a_offset + 0);
            *(b_offset + 1) = (i + 1 < m) ? *(a_offset + 1) : 0;
            *(b_offset + 2) = (i + 2 < m) ? *(a_offset + 2) : 0;
            *(b_offset + 3) = 0;
        }
    }
    return 0;
}
