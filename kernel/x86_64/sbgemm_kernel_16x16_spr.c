/***************************************************************************
 * Copyright (c) 2021, The OpenBLAS Project
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
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * *****************************************************************************/

#include <immintrin.h>
#include <string.h>
#include "common.h"

typedef struct {
	char palette_id;
	char start_row;
	char dummy0[14];  // bytes 2-15 reserved, must be zero
	short tile_colsb[8];
	char dummy1[16];  // bytes 32-47 reserved, must be zero
	char tile_rows[8];
	char dummy2[16];  // bytes 56-63 reserved, must be zero
} tilecfg;

/* tile0/tile1 -- A (m x 2k)
 * tile2/tile3 -- B (2k x n)
 * tile4-7 -- C (m x n)
 */
#define TCONF(cfg, m, n, k2) \
	memset(&cfg, 0, sizeof(tilecfg)); \
	cfg.palette_id = 1; \
	cfg.tile_rows[0] = m; \
	cfg.tile_rows[1] = m; \
	cfg.tile_rows[2] = k2>>1; \
	cfg.tile_rows[3] = k2>>1; \
	cfg.tile_rows[4] = m; \
	cfg.tile_rows[5] = m; \
	cfg.tile_rows[6] = m; \
	cfg.tile_rows[7] = m; \
	cfg.tile_colsb[0] = k2<<1; \
	cfg.tile_colsb[1] = k2<<1; \
	cfg.tile_colsb[2] = n * 4; \
	cfg.tile_colsb[3] = n * 4; \
	cfg.tile_colsb[4] = n * 4; \
	cfg.tile_colsb[5] = n * 4; \
	cfg.tile_colsb[6] = n * 4; \
	cfg.tile_colsb[7] = n * 4; \
	_tile_loadconfig(&cfg);

#define T_A0	0
#define T_A1	1
#define T_B0	2
#define T_B1	3
#define T_C00	4
#define T_C01	5
#define T_C10	6
#define T_C11	7

// FIXME: gcc11 seem have problem in tile load/store address calc,
// need to multiply with element size (2 or 4) here.
#define LOAD_A(M, N) _tile_loadd(T_A##M, ptr_a##M, lda * 2)
#define LOAD_A_TAIL(M, N) {\
	__m256i ymm = _mm256_loadu_epi16(ptr_a##M); \
	__m512i zmm = _mm512_cvtepu16_epi32(ymm); \
	_mm512_storeu_epi16(tail_a + 16 * M, zmm); \
	_tile_loadd(T_A##M, tail_a + 16 * 2 * M, 2 * 2); \
}
#define MASK_LOAD_A_TAIL(M, N) {\
	__m256i ymm = _mm256_maskz_loadu_epi16(amask, ptr_a##M); \
	__m512i zmm = _mm512_cvtepu16_epi32(ymm); \
	_mm512_storeu_epi16(tail_a + 16 * M, zmm); \
	_tile_loadd(T_A##M, tail_a + 16 * 2 * M, 2 * 2); \
}
#define LOAD_B(M, N) _tile_loadd(T_B##N, ptr_b##N, ldb * 2)
#define LOAD_B_TAIL(M, N) {\
	__m256i ymm = _mm256_loadu_epi16(ptr_b##N); \
	__m512i zmm = _mm512_cvtepu16_epi32(ymm); \
	_mm512_storeu_epi16(tail_b + 16 * N, zmm); \
	_tile_loadd(T_B##N, tail_b + 16 * 2 * N, 2 * 2); \
}
#define MASK_LOAD_B_TAIL(M, N) {\
	__m256i ymm = _mm256_maskz_loadu_epi16(bmask, ptr_b##N); \
	__m512i zmm = _mm512_cvtepu16_epi32(ymm); \
	_mm512_storeu_epi16(tail_b + 16 * N, zmm); \
	_tile_loadd(T_B##N, tail_b + 16 * 2 * N, 2 * 2); \
}

#define LOAD_C(M, N) _tile_loadd(T_C##M##N, ptr_c##M##N, ldc * 4)
#define MATMUL(M, N) _tile_dpbf16ps(T_C##M##N, T_A##M, T_B##N)
#define STORE_C(M, N) _tile_stored(T_C##M##N, ptr_c##M##N, ldc * 4)


int CNAME (BLASLONG im, BLASLONG in, BLASLONG k, FLOAT alpha, IFLOAT * iA, IFLOAT * iB, FLOAT * C, BLASLONG ldc)
{
	/* transport to Row Major matrix for AMX requirement */
	BLASLONG m, n;
	IFLOAT *A, *B;
	m = in;
	n = im;
	A = iB;
	B = iA;

	IFLOAT *ptr_a = A, *ptr_b = B;
	IFLOAT *ptr_b0, *ptr_b1;
	IFLOAT *ptr_a0, *ptr_a1;
	FLOAT *ptr_c = C;
	FLOAT *ptr_c00, *ptr_c01, *ptr_c10, *ptr_c11;

	BLASLONG lda, ldb;
	BLASLONG m_count = m;
	BLASLONG n_count, k_count;

	IFLOAT tail_a[32 * 2] __attribute__ ((aligned (64)));
	IFLOAT tail_b[32 * 2] __attribute__ ((aligned (64)));
	tilecfg cfg;

	for (; m_count > 31; m_count -= 32) {
		ptr_b = B;

		ptr_c00 = ptr_c;
		ptr_c01 = ptr_c00 + 16;
		ptr_c10 = ptr_c + 16 * ldc;
		ptr_c11 = ptr_c10 + 16;
		ptr_c += 32 * ldc;
		n_count = n;
		for (; n_count > 31; n_count -= 32) {
			ptr_a0 = ptr_a;
			ptr_a1 = ptr_a + 16 * k;

			ptr_b0 = ptr_b;
			ptr_b1 = ptr_b + 16 * k;
			ptr_b += 32 * k;

			lda = 32;
			ldb = 32;
			TCONF(cfg, 16, 16, 32);
			LOAD_C(0, 0); LOAD_C(0, 1);
			LOAD_C(1, 0); LOAD_C(1, 1);
			k_count = k;
			for (; k_count > 31; k_count -= 32) {
				LOAD_A(0, x); LOAD_A(1, x);
				ptr_a0 += 16 * 32;
				ptr_a1 += 16 * 32;
				LOAD_B(x, 0); LOAD_B(x, 1);
				ptr_b0 += 16 * 32;
				ptr_b1 += 16 * 32;

				MATMUL(0, 0); MATMUL(0, 1);
				MATMUL(1, 0); MATMUL(1, 1);
			}
			STORE_C(0, 0); STORE_C(0, 1);
			STORE_C(1, 0); STORE_C(1, 1);
			if (k_count > 1) {
				/* still have more than 2*k */
				int remain_k2 = k_count & ~1;
				k_count -= remain_k2;
				lda = remain_k2;
				TCONF(cfg, 16, 16, remain_k2);
				/* reconfig will clear all tiles,
				 * need to store/load again
				 */
				LOAD_C(0, 0); LOAD_C(0, 1);
				LOAD_C(1, 0); LOAD_C(1, 1);

				LOAD_A(0, x); LOAD_A(1, x);
				ptr_a0 += 16 * remain_k2;
				ptr_a1 += 16 * remain_k2;
				LOAD_B(x, 0); LOAD_B(x, 1);
				ptr_b0 += 16 * remain_k2;
				ptr_b1 += 16 * remain_k2;

				MATMUL(0, 0); MATMUL(0, 1);
				MATMUL(1, 0); MATMUL(1, 1);

				STORE_C(0, 0); STORE_C(0, 1);
				STORE_C(1, 0); STORE_C(1, 1);
			}
			if (k_count > 0) {
				/* still have odd tail k, need to transform into 2*k */
				TCONF(cfg, 16, 16, 2);

				LOAD_C(0, 0); LOAD_C(0, 1);
				LOAD_C(1, 0); LOAD_C(1, 1);

				LOAD_A_TAIL(0, x); LOAD_A_TAIL(1, x);
				LOAD_B_TAIL(x, 0); LOAD_B_TAIL(x, 1);

				MATMUL(0, 0); MATMUL(0, 1);
				MATMUL(1, 0); MATMUL(1, 1);

				STORE_C(0, 0); STORE_C(0, 1);
				STORE_C(1, 0); STORE_C(1, 1);
			}
			ptr_c00 += 32;
			ptr_c01 += 32;
			ptr_c10 += 32;
			ptr_c11 += 32;
		}
		for (; n_count > 0; n_count -= 16) {
			int tail_n = (n_count > 16) ? 16: n_count;
			__mmask16 bmask = (1UL << tail_n) - 1;
			ptr_a0 = ptr_a;
			ptr_a1 = ptr_a + 16 * k;

			ptr_b0 = ptr_b;
			ptr_b += tail_n * k;

			lda = 32;
			ldb = 2 * tail_n;
			TCONF(cfg, 16, tail_n, 32);
			LOAD_C(0, 0);
			LOAD_C(1, 0);
			k_count = k;
			for (; k_count > 31; k_count -= 32) {
				LOAD_A(0, x); LOAD_A(1, x);
				ptr_a0 += 16 * 32;
				ptr_a1 += 16 * 32;
				LOAD_B(x, 0);
				ptr_b0 += tail_n * 32;

				MATMUL(0, 0);
				MATMUL(1, 0);
			}
			STORE_C(0, 0);
			STORE_C(1, 0);
			if (k_count > 1) {
				/* still have more than 2*k */
				int remain_k2 = k_count & ~1;
				k_count -= remain_k2;
				lda = remain_k2;
				TCONF(cfg, 16, tail_n, remain_k2);
				/* reconfig will clear all tiles,
				 * need to store/load again
				 */
				LOAD_C(0, 0);
				LOAD_C(1, 0);

				LOAD_A(0, x); LOAD_A(1, x);
				ptr_a0 += 16 * remain_k2;
				ptr_a1 += 16 * remain_k2;
				LOAD_B(x, 0);
				ptr_b0 += tail_n * remain_k2;

				MATMUL(0, 0);
				MATMUL(1, 0);

				STORE_C(0, 0);
				STORE_C(1, 0);
			}
			if (k_count > 0) {
				/* still have odd tail k, need to transform into 2*k */
				TCONF(cfg, 16, tail_n, 2);

				LOAD_C(0, 0);
				LOAD_C(1, 0);

				LOAD_A_TAIL(0, x); LOAD_A_TAIL(1, x);
				MASK_LOAD_B_TAIL(x, 0);
				MATMUL(0, 0);
				MATMUL(1, 0);

				STORE_C(0, 0);
				STORE_C(1, 0);
			}
			ptr_c00 += tail_n;
			ptr_c10 += tail_n;
		}
		ptr_a += 32 * k;
	}
	for (; m_count > 0; m_count -= 16) {
		// process at most 16 m at a time
		int tail_m = (m_count > 16) ? 16: m_count;
		__mmask16 amask = (1UL << tail_m) - 1;

		ptr_b = B;

		ptr_c00 = ptr_c;
		ptr_c01 = ptr_c00 + 16;
		ptr_c += tail_m * ldc;
		n_count = n;
		for (; n_count > 31; n_count -= 32) {
			ptr_a0 = ptr_a;

			ptr_b0 = ptr_b;
			ptr_b1 = ptr_b + 16 * k;
			ptr_b += 32 * k;

			lda = 32;
			ldb = 32;
			TCONF(cfg, tail_m, 16, 32);
			LOAD_C(0, 0); LOAD_C(0, 1);
			k_count = k;
			for (; k_count > 31; k_count -= 32) {
				LOAD_A(0, x);
				ptr_a0 += tail_m * 32;
				LOAD_B(x, 0); LOAD_B(x, 1);
				ptr_b0 += 16 * 32;
				ptr_b1 += 16 * 32;

				MATMUL(0, 0); MATMUL(0, 1);
			}
			STORE_C(0, 0); STORE_C(0, 1);
			if (k_count > 1) {
				/* still have more than 2*k */
				int remain_k2 = k_count & ~1;
				k_count -= remain_k2;
				lda = remain_k2;
				TCONF(cfg, tail_m, 16, remain_k2);
				/* reconfig will clear all tiles,
				 * need to store/load again
				 */
				LOAD_C(0, 0); LOAD_C(0, 1);

				LOAD_A(0, x);
				ptr_a0 += tail_m * remain_k2;
				LOAD_B(x, 0); LOAD_B(x, 1);
				ptr_b0 += 16 * remain_k2;
				ptr_b1 += 16 * remain_k2;

				MATMUL(0, 0); MATMUL(0, 1);

				STORE_C(0, 0); STORE_C(0, 1);
			}
			if (k_count > 0) {
				/* still have odd tail k, need to transform into 2*k */
				TCONF(cfg, tail_m, 16, 2);

				LOAD_C(0, 0); LOAD_C(0, 1);

				MASK_LOAD_A_TAIL(0, x);
				LOAD_B_TAIL(x, 0); LOAD_B_TAIL(x, 1);

				MATMUL(0, 0); MATMUL(0, 1);

				STORE_C(0, 0); STORE_C(0, 1);
			}
			ptr_c00 += 32;
			ptr_c01 += 32;
		}
		for (; n_count > 0; n_count -= 16) {
			int tail_n = (n_count > 16) ? 16: n_count;
			__mmask16 bmask = (1UL << tail_n) - 1;
			ptr_a0 = ptr_a;

			ptr_b0 = ptr_b;
			ptr_b += tail_n * k;

			lda = 32;
			ldb = 2 * tail_n;
			TCONF(cfg, tail_m, tail_n, 32);
			LOAD_C(0, 0);
			k_count = k;
			for (; k_count > 31; k_count -= 32) {
				LOAD_A(0, x);
				ptr_a0 += tail_m * 32;
				LOAD_B(x, 0);
				ptr_b0 += tail_n * 32;

				MATMUL(0, 0);
			}
			STORE_C(0, 0);
			if (k_count > 1) {
				/* still have more than 2*k */
				int remain_k2 = k_count & ~1;
				k_count -= remain_k2;
				lda = remain_k2;
				TCONF(cfg, tail_m, tail_n, remain_k2);
				/* reconfig will clear all tiles,
				 * need to store/load again
				 */
				LOAD_C(0, 0);

				LOAD_A(0, x);
				ptr_a0 += tail_m * remain_k2;
				LOAD_B(x, 0);
				ptr_b0 += tail_n * remain_k2;

				MATMUL(0, 0);

				STORE_C(0, 0);
			}
			if (k_count > 0) {
				/* still have odd tail k, need to transform into 2*k */
				TCONF(cfg, tail_m, tail_n, 2);

				LOAD_C(0, 0);

				MASK_LOAD_A_TAIL(0, x);
				MASK_LOAD_B_TAIL(x, 0);
				MATMUL(0, 0);

				STORE_C(0, 0);
			}
			ptr_c00 += tail_n;
		}
		ptr_a += tail_m * k;
	}
	return 0;
}
