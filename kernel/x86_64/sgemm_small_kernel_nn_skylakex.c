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
#include <stdio.h>

#define DECLARE_RESULT_512(M, N) __m512 result##M##N = _mm512_setzero_ps()
#define LOAD_A_512(M, N) __m512 Aval##M = _mm512_loadu_ps(&A[lda * k + i + (M*16)])
#define MASK_LOAD_A_512(M, N) __m512 Aval##M = _mm512_maskz_loadu_ps(mask, &A[lda * k + i + (M*16)])
#define BROADCAST_LOAD_B_512(M, N) __m512 Bval##N = _mm512_broadcastss_ps(_mm_load_ss(&B[k + ldb * (j+N)]))
#define MATMUL_512(M, N) result##M##N = _mm512_fmadd_ps(Aval##M, Bval##N, result##M##N)
#if defined(B0)
#define STORE_512(M, N) result##M##N = _mm512_mul_ps(result##M##N, alpha_512); \
			_mm512_storeu_ps(&C[(j+N)*ldc + i + (M*16)], result##M##N)
#define MASK_STORE_512(M, N) result##M##N = _mm512_mul_ps(result##M##N, alpha_512); \
			_mm512_mask_storeu_ps(&C[(j+N)*ldc + i + (M*16)], mask, result##M##N)
#else
#define STORE_512(M, N) \
	BLASLONG offset##M##N = (j+N)*ldc + i + (M*16); \
	result##M##N = _mm512_mul_ps(result##M##N, alpha_512); \
	asm("vfmadd231ps (%1, %2, 4), %3, %0": "+v"(result##M##N):"r"(&C), "r"(offset##M##N), "v"(beta_512)); \
	_mm512_storeu_ps(&C[offset##M##N], result##M##N)
#define MASK_STORE_512(M, N) \
	BLASLONG offset##M##N = (j+N)*ldc + i + (M*16); \
	result##M##N = _mm512_mul_ps(result##M##N, alpha_512); \
	asm("vfmadd231ps (%1, %2, 4), %3, %0 %{%4%}": "+v"(result##M##N):"r"(&C), "r"(offset##M##N), "v"(beta_512), "k"(mask)); \
	_mm512_mask_storeu_ps(&C[offset##M##N], mask, result##M##N)
#endif

#define DECLARE_RESULT_256(M, N) __m256 result##M##N = _mm256_setzero_ps()
#define LOAD_A_256(M, N) __m256 Aval##M = _mm256_loadu_ps(&A[lda * k + i + (M*8)])
#define BROADCAST_LOAD_B_256(M, N) __m256 Bval##N = _mm256_broadcastss_ps(_mm_load_ss(&B[k + ldb * (j+N)]))
#define MATMUL_256(M, N) result##M##N = _mm256_fmadd_ps(Aval##M, Bval##N, result##M##N)
#if defined(B0)
#define STORE_256(M, N) result##M##N = _mm256_mul_ps(result##M##N, alpha_256); \
			_mm256_storeu_ps(&C[(j+N)*ldc + i + (M*8)], result##M##N)
#else
#define STORE_256(M, N) \
	BLASLONG offset##M##N = (j+N)*ldc + i + (M*8); \
	result##M##N = _mm256_mul_ps(result##M##N, alpha_256); \
	asm("vfmadd231ps (%1, %2, 4), %3, %0": "+v"(result##M##N):"r"(&C), "r"(offset##M##N), "v"(beta_256)); \
	_mm256_storeu_ps(&C[offset##M##N], result##M##N)
#endif

#define DECLARE_RESULT_128(M, N) __m128 result##M##N; asm("vpxorq %0, %0, %0": "+v"(result##M##N):)
#define LOAD_A_128(M, N) __m128 Aval##M = _mm_maskz_loadu_ps(mask, &A[lda * k + i + (M*4)])
#define BROADCAST_LOAD_B_128(M, N) __m128 Bval##N = _mm_broadcastss_ps(_mm_load_ss(&B[k + ldb * (j+N)]))
#define MATMUL_128(M, N) result##M##N = _mm_fmadd_ps(Aval##M, Bval##N, result##M##N)
#if defined(B0)
#define STORE_128(M, N) result##M##N = _mm_maskz_mul_ps(mask, result##M##N, alpha_128); \
			_mm_mask_storeu_ps(&C[(j+N)*ldc + i + (M*4)], mask, result##M##N)
#else
#define STORE_128(M, N) \
	BLASLONG offset##M##N = (j+N)*ldc + i + (M*4); \
	result##M##N = _mm_maskz_mul_ps(mask, result##M##N, alpha_128); \
	asm("vfmadd231ps (%1, %2, 4), %3, %0": "+v"(result##M##N):"r"(&C), "r"(offset##M##N), "v"(beta_128)); \
	_mm_mask_storeu_ps(&C[offset##M##N], mask, result##M##N)
#endif

#define DECLARE_RESULT_S(M, N) float result##M##N = 0;
#define LOAD_A_S(M, N) float Aval##M = A[lda * k + i + M]
#define BROADCAST_LOAD_B_S(M, N) float Bval##N = B[k + ldb * (j+N)]
#define MATMUL_S(M, N) result##M##N += Aval##M * Bval##N
#if defined(B0)
#define STORE_S(M, N) C[(j+N)*ldc + i + M] = result##M##N * alpha
#else
#define STORE_S(M, N) C[(j+N)*ldc + i + M] = result##M##N * alpha + C[(j+N)*ldc + i + M] * beta
#endif

#if defined(B0)
int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT * C, BLASLONG ldc)
#else
int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT beta, FLOAT * C, BLASLONG ldc)
#endif
{
	// column major
	BLASLONG i, j, k;

	BLASLONG m64 = M & ~63;
	BLASLONG m32 = M & ~31;
	BLASLONG m16 = M & ~15;
	BLASLONG m8 = M & ~7;
	BLASLONG m4 = M & ~3;
	BLASLONG m2 = M & ~1;

	BLASLONG n4 = N & ~3;
	BLASLONG n2 = N & ~1;

	__mmask8 mask = 0xff;  // just use to avoid SSE instruction

	__m512 alpha_512 = _mm512_broadcastss_ps(_mm_load_ss(&alpha));
#if !defined(B0)
	__m512 beta_512 = _mm512_broadcastss_ps(_mm_load_ss(&beta));
#endif

	for (i = 0; i < m64; i += 64) {
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
			DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1); DECLARE_RESULT_512(2, 1); DECLARE_RESULT_512(3, 1);
			DECLARE_RESULT_512(0, 2); DECLARE_RESULT_512(1, 2); DECLARE_RESULT_512(2, 2); DECLARE_RESULT_512(3, 2);
			DECLARE_RESULT_512(0, 3); DECLARE_RESULT_512(1, 3); DECLARE_RESULT_512(2, 3); DECLARE_RESULT_512(3, 3);

			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x); LOAD_A_512(1, x); LOAD_A_512(2, x); LOAD_A_512(3, x);

				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				BROADCAST_LOAD_B_512(x, 2); BROADCAST_LOAD_B_512(x, 3);

				MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
				MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
				MATMUL_512(0, 2); MATMUL_512(1, 2); MATMUL_512(2, 2); MATMUL_512(3, 2);
				MATMUL_512(0, 3); MATMUL_512(1, 3); MATMUL_512(2, 3); MATMUL_512(3, 3);
			}
			STORE_512(0, 0); STORE_512(1, 0); STORE_512(2, 0); STORE_512(3, 0);
			STORE_512(0, 1); STORE_512(1, 1); STORE_512(2, 1); STORE_512(3, 1);
			STORE_512(0, 2); STORE_512(1, 2); STORE_512(2, 2); STORE_512(3, 2);
			STORE_512(0, 3); STORE_512(1, 3); STORE_512(2, 3); STORE_512(3, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
			DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1); DECLARE_RESULT_512(2, 1); DECLARE_RESULT_512(3, 1);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x); LOAD_A_512(1, x); LOAD_A_512(2, x); LOAD_A_512(3, x);
				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
				MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
			}
			STORE_512(0, 0); STORE_512(1, 0); STORE_512(2, 0); STORE_512(3, 0);
			STORE_512(0, 1); STORE_512(1, 1); STORE_512(2, 1); STORE_512(3, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x); LOAD_A_512(1, x); LOAD_A_512(2, x); LOAD_A_512(3, x);
				BROADCAST_LOAD_B_512(x, 0);
				MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
			}
			STORE_512(0, 0); STORE_512(1, 0); STORE_512(2, 0); STORE_512(3, 0);
		}
	}
	for (; i < m32; i += 32) {
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0);
			DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1);
			DECLARE_RESULT_512(0, 2); DECLARE_RESULT_512(1, 2);
			DECLARE_RESULT_512(0, 3); DECLARE_RESULT_512(1, 3);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x); LOAD_A_512(1, x);
				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				BROADCAST_LOAD_B_512(x, 2); BROADCAST_LOAD_B_512(x, 3);

				MATMUL_512(0, 0); MATMUL_512(1, 0);
				MATMUL_512(0, 1); MATMUL_512(1, 1);
				MATMUL_512(0, 2); MATMUL_512(1, 2);
				MATMUL_512(0, 3); MATMUL_512(1, 3);
			}
			STORE_512(0, 0); STORE_512(1, 0);
			STORE_512(0, 1); STORE_512(1, 1);
			STORE_512(0, 2); STORE_512(1, 2);
			STORE_512(0, 3); STORE_512(1, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0);
			DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x); LOAD_A_512(1, x);
				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				MATMUL_512(0, 0); MATMUL_512(1, 0);
				MATMUL_512(0, 1); MATMUL_512(1, 1);
			}
			STORE_512(0, 0); STORE_512(1, 0);
			STORE_512(0, 1); STORE_512(1, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x); LOAD_A_512(1, x);
				BROADCAST_LOAD_B_512(x, 0);
				MATMUL_512(0, 0); MATMUL_512(1, 0);
			}
			STORE_512(0, 0); STORE_512(1, 0);
		}
	}
	for (; i < m16; i += 16) {
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_512(0, 0);
			DECLARE_RESULT_512(0, 1);
			DECLARE_RESULT_512(0, 2);
			DECLARE_RESULT_512(0, 3);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x);
				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				BROADCAST_LOAD_B_512(x, 2); BROADCAST_LOAD_B_512(x, 3);

				MATMUL_512(0, 0);
				MATMUL_512(0, 1);
				MATMUL_512(0, 2);
				MATMUL_512(0, 3);
			}
			STORE_512(0, 0);
			STORE_512(0, 1);
			STORE_512(0, 2);
			STORE_512(0, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_512(0, 0);
			DECLARE_RESULT_512(0, 1);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x);
				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				MATMUL_512(0, 0);
				MATMUL_512(0, 1);
			}
			STORE_512(0, 0);
			STORE_512(0, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_512(0, 0);
			for (k = 0; k < K; k++) {
				LOAD_A_512(0, x);
				BROADCAST_LOAD_B_512(x, 0);
				MATMUL_512(0, 0);
			}
			STORE_512(0, 0);
		}
	}
	if (M - i > 0) {
		register __mmask16 mask asm("k1") = (1UL << (M - i)) - 1;
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_512(0, 0);
			DECLARE_RESULT_512(0, 1);
			DECLARE_RESULT_512(0, 2);
			DECLARE_RESULT_512(0, 3);
			for (k = 0; k < K; k++) {
				MASK_LOAD_A_512(0, x);
				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				BROADCAST_LOAD_B_512(x, 2); BROADCAST_LOAD_B_512(x, 3);

				MATMUL_512(0, 0);
				MATMUL_512(0, 1);
				MATMUL_512(0, 2);
				MATMUL_512(0, 3);
			}
			MASK_STORE_512(0, 0);
			MASK_STORE_512(0, 1);
			MASK_STORE_512(0, 2);
			MASK_STORE_512(0, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_512(0, 0);
			DECLARE_RESULT_512(0, 1);
			for (k = 0; k < K; k++) {
				MASK_LOAD_A_512(0, x);
				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
				MATMUL_512(0, 0);
				MATMUL_512(0, 1);
			}
			MASK_STORE_512(0, 0);
			MASK_STORE_512(0, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_512(0, 0);
			for (k = 0; k < K; k++) {
				MASK_LOAD_A_512(0, x);
				BROADCAST_LOAD_B_512(x, 0);
				MATMUL_512(0, 0);
			}
			MASK_STORE_512(0, 0);
		}
		return;
	}
	__m256 alpha_256 = _mm256_broadcastss_ps(_mm_load_ss(&alpha));
#if !defined(B0)
	__m256 beta_256 = _mm256_broadcastss_ps(_mm_load_ss(&beta));
#endif
	for (; i < m8; i += 8) {
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_256(0, 0);
			DECLARE_RESULT_256(0, 1);
			DECLARE_RESULT_256(0, 2);
			DECLARE_RESULT_256(0, 3);
			for (k = 0; k < K; k++) {
				LOAD_A_256(0, x);
				BROADCAST_LOAD_B_256(x, 0); BROADCAST_LOAD_B_256(x, 1);
				BROADCAST_LOAD_B_256(x, 2); BROADCAST_LOAD_B_256(x, 3);

				MATMUL_256(0, 0);
				MATMUL_256(0, 1);
				MATMUL_256(0, 2);
				MATMUL_256(0, 3);
			}
			STORE_256(0, 0);
			STORE_256(0, 1);
			STORE_256(0, 2);
			STORE_256(0, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_256(0, 0);
			DECLARE_RESULT_256(0, 1);
			for (k = 0; k < K; k++) {
				LOAD_A_256(0, x);
				BROADCAST_LOAD_B_256(x, 0); BROADCAST_LOAD_B_256(x, 1);
				MATMUL_256(0, 0);
				MATMUL_256(0, 1);
			}
			STORE_256(0, 0);
			STORE_256(0, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_256(0, 0);
			for (k = 0; k < K; k++) {
				LOAD_A_256(0, x);
				BROADCAST_LOAD_B_256(x, 0);
				MATMUL_256(0, 0);
			}
			STORE_256(0, 0);
		}
	}
	__m128 alpha_128 = _mm_broadcastss_ps(_mm_load_ss(&alpha));
#if !defined(B0)
	__m128 beta_128 = _mm_broadcastss_ps(_mm_load_ss(&beta));
#endif
	for (; i < m4; i += 4) {
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_128(0, 0);
			DECLARE_RESULT_128(0, 1);
			DECLARE_RESULT_128(0, 2);
			DECLARE_RESULT_128(0, 3);
			for (k = 0; k < K; k++) {
				LOAD_A_128(0, x);
				BROADCAST_LOAD_B_128(x, 0); BROADCAST_LOAD_B_128(x, 1);
				BROADCAST_LOAD_B_128(x, 2); BROADCAST_LOAD_B_128(x, 3);

				MATMUL_128(0, 0);
				MATMUL_128(0, 1);
				MATMUL_128(0, 2);
				MATMUL_128(0, 3);
			}
			STORE_128(0, 0);
			STORE_128(0, 1);
			STORE_128(0, 2);
			STORE_128(0, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_128(0, 0);
			DECLARE_RESULT_128(0, 1);
			for (k = 0; k < K; k++) {
				LOAD_A_128(0, x);
				BROADCAST_LOAD_B_128(x, 0); BROADCAST_LOAD_B_128(x, 1);
				MATMUL_128(0, 0);
				MATMUL_128(0, 1);
			}
			STORE_128(0, 0);
			STORE_128(0, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_128(0, 0);
			for (k = 0; k < K; k++) {
				LOAD_A_128(0, x);
				BROADCAST_LOAD_B_128(x, 0);
				MATMUL_128(0, 0);
			}
			STORE_128(0, 0);
		}
	}
	for (; i < m2; i += 2) {
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_S(0, 0); DECLARE_RESULT_S(1, 0);
			DECLARE_RESULT_S(0, 1); DECLARE_RESULT_S(1, 1);
			DECLARE_RESULT_S(0, 2); DECLARE_RESULT_S(1, 2);
			DECLARE_RESULT_S(0, 3); DECLARE_RESULT_S(1, 3);
			for (k = 0; k < K; k++) {
				LOAD_A_S(0, x); LOAD_A_S(1, x);
				BROADCAST_LOAD_B_S(x, 0); BROADCAST_LOAD_B_S(x, 1);
				BROADCAST_LOAD_B_S(x, 2); BROADCAST_LOAD_B_S(x, 3);

				MATMUL_S(0, 0); MATMUL_S(1, 0);
				MATMUL_S(0, 1); MATMUL_S(1, 1);
				MATMUL_S(0, 2); MATMUL_S(1, 2);
				MATMUL_S(0, 3); MATMUL_S(1, 3);
			}
			STORE_S(0, 0); STORE_S(1, 0);
			STORE_S(0, 1); STORE_S(1, 1);
			STORE_S(0, 2); STORE_S(1, 2);
			STORE_S(0, 3); STORE_S(1, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_S(0, 0); DECLARE_RESULT_S(1, 0);
			DECLARE_RESULT_S(0, 1); DECLARE_RESULT_S(1, 1);
			for (k = 0; k < K; k++) {
				LOAD_A_S(0, x); LOAD_A_S(1, x);
				BROADCAST_LOAD_B_S(x, 0); BROADCAST_LOAD_B_S(x, 1);
				MATMUL_S(0, 0); MATMUL_S(1, 0);
				MATMUL_S(0, 1); MATMUL_S(1, 1);
			}
			STORE_S(0, 0); STORE_S(1, 0);
			STORE_S(0, 1); STORE_S(1, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_S(0, 0); DECLARE_RESULT_S(1, 0);
			for (k = 0; k < K; k++) {
				LOAD_A_S(0, x); LOAD_A_S(1, x);
				BROADCAST_LOAD_B_S(x, 0);
				MATMUL_S(0, 0); MATMUL_S(1, 0);
			}
			STORE_S(0, 0); STORE_S(1, 0);
		}
	}
	for (; i < M; i += 1) {
		for (j = 0; j < n4; j += 4) {
			DECLARE_RESULT_S(0, 0);
			DECLARE_RESULT_S(0, 1);
			DECLARE_RESULT_S(0, 2);
			DECLARE_RESULT_S(0, 3);
			for (k = 0; k < K; k++) {
				LOAD_A_S(0, x);
				BROADCAST_LOAD_B_S(x, 0); BROADCAST_LOAD_B_S(x, 1);
				BROADCAST_LOAD_B_S(x, 2); BROADCAST_LOAD_B_S(x, 3);

				MATMUL_S(0, 0);
				MATMUL_S(0, 1);
				MATMUL_S(0, 2);
				MATMUL_S(0, 3);
			}
			STORE_S(0, 0);
			STORE_S(0, 1);
			STORE_S(0, 2);
			STORE_S(0, 3);
		}
		for (; j < n2; j += 2) {
			DECLARE_RESULT_S(0, 0);
			DECLARE_RESULT_S(0, 1);
			for (k = 0; k < K; k++) {
				LOAD_A_S(0, x);
				BROADCAST_LOAD_B_S(x, 0); BROADCAST_LOAD_B_S(x, 1);
				MATMUL_S(0, 0);
				MATMUL_S(0, 1);
			}
			STORE_S(0, 0);
			STORE_S(0, 1);
		}
		for (; j < N; j++) {
			DECLARE_RESULT_S(0, 0);
			for (k = 0; k < K; k++) {
				LOAD_A_S(0, x); LOAD_A_S(1, x);
				BROADCAST_LOAD_B_S(x, 0);
				MATMUL_S(0, 0);
			}
			STORE_S(0, 0);
		}
	}
}
