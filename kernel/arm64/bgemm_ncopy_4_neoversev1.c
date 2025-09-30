/***************************************************************************
 * Copyright (c) 2025, The OpenBLAS Project
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
#include <arm_neon.h>

#include "common.h"

int CNAME(BLASLONG n, BLASLONG m, IFLOAT *input, BLASLONG lda, IFLOAT *output) {
    const int num_accumulators = 4;
    const int m_accumulators = m & -4;

    const int n4 = n & -4;
    const int n_rest = n - n4;

    const int m_rest = m - m_accumulators;

    for (size_t m_step = 0; m_step < m_accumulators; m_step += num_accumulators) {
        const uint16_t* inner_input = input;

        // Potential for vld1q here with transpose
        for (size_t n_step = 0; n_step < n4; n_step += 4) {
            uint16x4_t a_vec0 = vld1_u16(inner_input + 0 * lda);
            uint16x4_t a_vec1 = vld1_u16(inner_input + 1 * lda);
            uint16x4_t a_vec2 = vld1_u16(inner_input + 2 * lda);
            uint16x4_t a_vec3 = vld1_u16(inner_input + 3 * lda);

            vst1_u16(output, a_vec0);
            vst1_u16(output + 4, a_vec1);
            vst1_u16(output + 8, a_vec2);
            vst1_u16(output + 12, a_vec3);

            output += 16;
            inner_input += 4;
        }

        if (n_rest) {
            for (BLASLONG line = 0; line < num_accumulators; line++) {
                output[0] = inner_input[0];
                output[1] = n_rest == 1 ? 0 : inner_input[1];
                output[2] = n_rest <= 2 ? 0 : inner_input[2];
                output[3] = n_rest <= 3 ? 0 : inner_input[3];

                inner_input += lda;
                output += 4;
            }
        }

        input += lda * num_accumulators;
    }

    if (m_rest & 2) {
        const uint16_t* inner_input = input;
        for (size_t n_step = 0; n_step < n4; n_step += 4) {
            uint16x4_t a_vec0 = vld1_u16(inner_input);
            uint16x4_t a_vec1 = vld1_u16(inner_input + lda);

            vst1_u16(output, a_vec0);
            vst1_u16(output + 4, a_vec1);

            inner_input += 4;
            output += 8;
        }

        if (n_rest) {
            for (BLASLONG line = 0; line < 2; line++) {
                output[0] = inner_input[0];
                output[1] = n_rest == 1 ? 0 : inner_input[1];
                output[2] = n_rest <= 2 ? 0 : inner_input[2];
                output[3] = n_rest <= 3 ? 0 : inner_input[3];

                inner_input += lda;
                output += 4;
            }
        }

        input += lda * 2;
    }

    if (m_rest & 1) {
        for (size_t n_step = 0; n_step < n4; n_step += 4) {
            uint16x4_t a_vec0 = vld1_u16(input);

            vst1_u16(output, a_vec0);

            input += 4;
            output += 4;
        }

        if (n_rest) {
            output[0] = input[0];
            output[1] = n_rest == 1 ? 0 : input[1];
            output[2] = n_rest <= 2 ? 0 : input[2];
            output[3] = n_rest <= 3 ? 0 : input[3];
        }
    }

    return 0;
}
