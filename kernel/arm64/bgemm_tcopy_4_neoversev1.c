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

int CNAME(BLASLONG m, BLASLONG n, IFLOAT *input, BLASLONG lda, IFLOAT *output) {
    const int num_accumulators = 4;
    const int n_accumulators = (n & -num_accumulators);
    const int n_rest = n - n_accumulators;

    const int m4 = m & -4;
    const int m_rest = m - m4;

    for (size_t n_step = 0; n_step < n_accumulators; n_step += num_accumulators) {
        const uint16_t* inner_input = input;

        // Full 4x4 item transposes down the M dimension
        for (size_t m_step = 0; m_step < m4; m_step += 4) {
            // Load 4x4 block
            uint16x4_t a_vec0 = vld1_u16(inner_input);
            uint16x4_t a_vec1 = vld1_u16(inner_input + lda);
            uint16x4_t a_vec2 = vld1_u16(inner_input + 2 * lda);
            uint16x4_t a_vec3 = vld1_u16(inner_input + 3 * lda);

            // Transpose 4x4 blocks
            uint16x4_t out_vec0 = vzip1_u16(a_vec0, a_vec1);
            uint16x4_t out_vec1 = vzip2_u16(a_vec0, a_vec1);
            uint16x4_t out_vec2 = vzip1_u16(a_vec2, a_vec3);
            uint16x4_t out_vec3 = vzip2_u16(a_vec2, a_vec3);

            // Transpose 8x4 blocks
            a_vec0 = vreinterpret_u16_u32(vzip1_u32(vreinterpret_u32_u16(out_vec0), vreinterpret_u32_u16(out_vec2)));
            a_vec1 = vreinterpret_u16_u32(vzip2_u32(vreinterpret_u32_u16(out_vec0), vreinterpret_u32_u16(out_vec2)));
            a_vec2 = vreinterpret_u16_u32(vzip1_u32(vreinterpret_u32_u16(out_vec1), vreinterpret_u32_u16(out_vec3)));
            a_vec3 = vreinterpret_u16_u32(vzip2_u32(vreinterpret_u32_u16(out_vec1), vreinterpret_u32_u16(out_vec3)));

            vst1_u16(output, a_vec0);
            vst1_u16(output + 4, a_vec1);
            vst1_u16(output + 8, a_vec2);
            vst1_u16(output + 12, a_vec3);

            inner_input += 4 * lda;
            output += 16;
        }

        if (m_rest) {
            for (BLASLONG line = 0; line < num_accumulators; line++) {
                output[0] = inner_input[0];
                output[1] = m_rest == 1 ? 0 : *(inner_input + lda);
                output[2] = m_rest <= 2 ? 0 : *(inner_input + 2 * lda);
                output[3] = m_rest <= 3 ? 0 : *(inner_input + 3 * lda);

                inner_input++;
                output += 4;
            }
        }

        input += num_accumulators;
    }

    // Extract two remaining rows as 128-bit vector paired
    if (n_rest & 2) {
        const uint16_t* inner_input = input;
        for (size_t m_step = 0; m_step < m4; m_step += 4) {
            for (BLASLONG line = 0; line < 2; line++) {
                output[0] = *(inner_input + line);
                output[1] = *(inner_input + line + lda);
                output[2] = *(inner_input + line + 2 * lda);
                output[3] = *(inner_input + line + 3 * lda);

                output += 4;
            }

            inner_input += 4 * lda;
        }

        if (m_rest) {
            for (BLASLONG line = 0; line < 2; line++) {
                output[0] = *(inner_input + line);
                output[1] = m_rest == 1 ? 0 : *(inner_input + line + lda);
                output[2] = m_rest <= 2 ? 0 : *(inner_input + line + 2 * lda);
                output[3] = m_rest <= 3 ? 0 : *(inner_input + line + 3 * lda);

                output += 4;
            }
        }

        input += 2;
    }

    // Flatten final row
    if (n_rest & 1) {
        const uint16_t* inner_input = input;
        for (size_t m_step = 0; m_step < m4; m_step += 4) {
            output[0] = *inner_input;
            output[1] = *(inner_input + lda);
            output[2] = *(inner_input + 2 * lda);
            output[3] = *(inner_input + 3 * lda);

            inner_input += 4 * lda;
            output += 4;
        }

        if (m_rest) {
            output[0] = inner_input[0];
            output[1] = m_rest == 1 ? 0 : *(inner_input + lda);
            output[2] = m_rest <= 2 ? 0 : *(inner_input + 2 * lda);
            output[3] = m_rest <= 3 ? 0 : *(inner_input + 3 * lda);

            output += 4;
        }
    }

    return 0;
}
