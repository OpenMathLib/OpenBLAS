/***************************************************************************
Copyright (c) 2021, The OpenBLAS Project
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

#include <immintrin.h>
#include "common.h"

#define VMOVLDUP(addr, zmm) asm("vmovsldup (%1), %0": "=v"(zmm): "r"(addr))
#define VMOVHDUP(addr, zmm) asm("vmovshdup (%1), %0": "=v"(zmm): "r"(addr))
#define BROADCAST64(base, step, n, offset, zmm) \
	if (n == 0) asm("vbroadcastsd %2(%1), %0": "=v"(zmm): "r"(base), "n"(offset*2)); \
	else asm("vbroadcastsd %4(%1, %2, %3), %0": "=v"(zmm): "r"(base), "r"(step), "n"(n*2), "n"(offset*2))

#define DECLARE_A_PAIR(A) \
	__m512i A_lo_##A; __m512i A_hi_##A;

#define LOAD_A_PAIR(A) \
	VMOVLDUP(ptr_a##A, A_lo_##A); \
	VMOVHDUP(ptr_a##A, A_hi_##A);

#define LOAD_A_PAIR_TAIL(A) { \
	__m256i ymm = _mm256_loadu_si256(ptr_a##A); \
	__m512 zmm = (__m512) _mm512_cvtepu16_epi32(ymm); \
	A_lo_##A = (__m512i) _mm512_moveldup_ps(zmm); \
	A_hi_##A = (__m512i) _mm512_movehdup_ps(zmm); \
}

#define DECLARE_B_PAIR() \
	__m512i B_lo; __m512i B_hi;

#define BROADCAST_B_PAIR(Bx, By) \
	BROADCAST64(ptr_b##Bx, n_blksize, By, 0, B_lo); \
	BROADCAST64(ptr_b##Bx, n_blksize, By, 2, B_hi);

#define BROADCAST_B_PAIR_TAIL(Bx, By) {\
	__m128i xmm = (__m128i) _mm_load_sd(ptr_b##Bx + n_blksize * By); \
	xmm = _mm_cvtepu16_epi32(xmm); \
	B_lo = _mm512_broadcastd_epi32(xmm); \
	B_hi = _mm512_broadcastd_epi32((__m128i) _mm_permute_pd((__m128d) xmm, 0x1)); \
}

#define DECLARE_RESULT_4X(A, Bx, By) \
	__m512 result_00_##A##Bx##By = _mm512_setzero_ps(); \
	__m512 result_01_##A##Bx##By = _mm512_setzero_ps(); \
	__m512 result_10_##A##Bx##By = _mm512_setzero_ps(); \
	__m512 result_11_##A##Bx##By = _mm512_setzero_ps();

#define FMA(a, b, r) r = _mm512_dpbf16_ps(r, (__m512bh)a, (__m512bh)b)

#define MATMUL_4X(A, Bx, By) \
	FMA(A_lo_##A, B_lo, result_00_##A##Bx##By); \
	FMA(A_hi_##A, B_lo, result_01_##A##Bx##By); \
	FMA(A_lo_##A, B_hi, result_10_##A##Bx##By); \
	FMA(A_hi_##A, B_hi, result_11_##A##Bx##By);

#define STORE_4X(A, Bx, By)



int CNAME (BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, IFLOAT * A, IFLOAT * B, FLOAT * C, BLASLONG ldc)
{
	IFLOAT *ptr_a = A, *ptr_b = B, *ptr_c = C;
	IFLOAT *ptr_b0, *ptr_b1;
	IFLOAT *ptr_a0, *ptr_a1;
	BLASLONG n_count = n;
	BLASLONG m_count, k_count;
	BLASLONG n_blksize = 4 * k;

	for (; n_count > 23; n_count -= 24) {
		m_count = m;
		ptr_b0 = ptr_b;
		ptr_b1 = ptr_b0 + n_blksize * 3;
		for (; m_count > 15; m_count -= 16) {
			DECLARE_A_PAIR(0); DECLARE_B_PAIR();
			DECLARE_RESULT_4X(0, 0, 0); DECLARE_RESULT_4X(0, 0, 1); DECLARE_RESULT_4X(0, 0, 2);
			DECLARE_RESULT_4X(0, 1, 0); DECLARE_RESULT_4X(0, 1, 1); DECLARE_RESULT_4X(0, 1, 2);
			for (k_count = k; k_count > 1; k_count -=2) {
				LOAD_A_PAIR(0);
				BROADCAST_B_PAIR(0, 0); MATMUL_4X(0, 0, 0);
				BROADCAST_B_PAIR(0, 1); MATMUL_4X(0, 0, 1);
				BROADCAST_B_PAIR(0, 2); MATMUL_4X(0, 0, 2);
				BROADCAST_B_PAIR(1, 0); MATMUL_4X(0, 1, 0);
				BROADCAST_B_PAIR(1, 1); MATMUL_4X(0, 1, 1);
				BROADCAST_B_PAIR(1, 2); MATMUL_4X(0, 1, 2);
				ptr_b0 += 24 * 2;
				ptr_b1 += 24 * 2;
				ptr_a0 += 16 * 2;
			}
			if (k_count > 0) {
				LOAD_A_PAIR_TAIL(0);
				BROADCAST_B_PAIR_TAIL(0, 0); MATMUL_4X(0, 0, 0);
				BROADCAST_B_PAIR_TAIL(0, 1); MATMUL_4X(0, 0, 1);
				BROADCAST_B_PAIR_TAIL(0, 2); MATMUL_4X(0, 0, 2);
				BROADCAST_B_PAIR_TAIL(1, 0); MATMUL_4X(0, 1, 0);
				BROADCAST_B_PAIR_TAIL(1, 1); MATMUL_4X(0, 1, 1);
				BROADCAST_B_PAIR_TAIL(1, 2); MATMUL_4X(0, 1, 2);
				ptr_b0 += 24;
				ptr_b1 += 24;
				ptr_a0 += 16;
			}
		}
	}
}
