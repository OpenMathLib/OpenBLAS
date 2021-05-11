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
#include <memory.h>

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

#define LOAD_KA_512(M, N) __m512 Aval##M = _mm512_loadu_ps(&mbuf[(mi + M)*K + k]);
#define LOAD_KB_512(M, N) __m512 Bval##N = _mm512_loadu_ps(&B[(j + N)*ldb + k])
#define MASK_LOAD_KA_512(M, N) __m512 Aval##M = _mm512_maskz_loadu_ps(mask, &mbuf[(mi + M)*K + k])
#define MASK_LOAD_KB_512(M, N) __m512 Bval##N = _mm512_maskz_loadu_ps(mask, &B[(j + N)*ldb + k])
#if defined(B0)
#define STORE_REDUCE(M, N) C[(j+N)*ldc + i + M] = alpha * _mm512_reduce_add_ps(result##M##N);
#else
#define STORE_REDUCE(M, N) C[(j+N)*ldc + i + M] = alpha * _mm512_reduce_add_ps(result##M##N) + beta * C[(j+N)*ldc + i + M];
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
	if (M - i > 8) {
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
	int mm = M - i;
	if (mm) {
		FLOAT *mbuf = (FLOAT *) malloc(sizeof(FLOAT)*mm*K);
		__mmask8 mask8 = (1UL << mm) - 1;
		__mmask16 mask;
		BLASLONG k16 = K & ~15;
		BLASLONG k8 = K & ~7;
		for (k = 0; k < k8; k += 8) {
			__m256  r0, r1, r2, r3, r4, r5, r6, r7;
			__m256  t0, t1, t2, t3, t4, t5, t6, t7;
			r0 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(0 + k)]);
			r1 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(1 + k)]);
			r2 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(2 + k)]);
			r3 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(3 + k)]);
			r4 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(4 + k)]);
			r5 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(5 + k)]);
			r6 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(6 + k)]);
			r7 = _mm256_maskz_loadu_ps(mask8, &A[i + lda*(7 + k)]);

			t0 = _mm256_unpacklo_ps(r0, r1);
			t1 = _mm256_unpackhi_ps(r0, r1);
			t2 = _mm256_unpacklo_ps(r2, r3);
			t3 = _mm256_unpackhi_ps(r2, r3);
			t4 = _mm256_unpacklo_ps(r4, r5);
			t5 = _mm256_unpackhi_ps(r4, r5);
			t6 = _mm256_unpacklo_ps(r6, r7);
			t7 = _mm256_unpackhi_ps(r6, r7);

			r0 = _mm256_shuffle_ps(t0,t2,_MM_SHUFFLE(1,0,1,0));
			r1 = _mm256_shuffle_ps(t0,t2,_MM_SHUFFLE(3,2,3,2));
			r2 = _mm256_shuffle_ps(t1,t3,_MM_SHUFFLE(1,0,1,0));
			r3 = _mm256_shuffle_ps(t1,t3,_MM_SHUFFLE(3,2,3,2));
			r4 = _mm256_shuffle_ps(t4,t6,_MM_SHUFFLE(1,0,1,0));
			r5 = _mm256_shuffle_ps(t4,t6,_MM_SHUFFLE(3,2,3,2));
			r6 = _mm256_shuffle_ps(t5,t7,_MM_SHUFFLE(1,0,1,0));
			r7 = _mm256_shuffle_ps(t5,t7,_MM_SHUFFLE(3,2,3,2));

			t0 = _mm256_permute2f128_ps(r0, r4, 0x20);
			t1 = _mm256_permute2f128_ps(r1, r5, 0x20);
			t2 = _mm256_permute2f128_ps(r2, r6, 0x20);
			t3 = _mm256_permute2f128_ps(r3, r7, 0x20);
			t4 = _mm256_permute2f128_ps(r0, r4, 0x31);
			t5 = _mm256_permute2f128_ps(r1, r5, 0x31);
			t6 = _mm256_permute2f128_ps(r2, r6, 0x31);
			t7 = _mm256_permute2f128_ps(r3, r7, 0x31);

			switch (mm) {
				case 8: _mm256_storeu_ps(&mbuf[k + 7*K], t7);
				case 7: _mm256_storeu_ps(&mbuf[k + 6*K], t6);
				case 6: _mm256_storeu_ps(&mbuf[k + 5*K], t5);
				case 5: _mm256_storeu_ps(&mbuf[k + 4*K], t4);
				case 4: _mm256_storeu_ps(&mbuf[k + 3*K], t3);
				case 3: _mm256_storeu_ps(&mbuf[k + 2*K], t2);
				case 2: _mm256_storeu_ps(&mbuf[k + 1*K], t1);
				case 1: _mm256_storeu_ps(&mbuf[k + 0*K], t0);
			}
		}
		for (; k < K; k++) {
			for (int ii = 0; ii < mm; ii++) {
				mbuf[k + ii*K] = A[i + lda*k + ii];
			}
		}
		int mi = 0;
		for (; i < m4; i += 4, mi += 4) {
			for (j = 0; j < n4; j += 4) {
				DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
				DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1); DECLARE_RESULT_512(2, 1); DECLARE_RESULT_512(3, 1);
				DECLARE_RESULT_512(0, 2); DECLARE_RESULT_512(1, 2); DECLARE_RESULT_512(2, 2); DECLARE_RESULT_512(3, 2);
				DECLARE_RESULT_512(0, 3); DECLARE_RESULT_512(1, 3); DECLARE_RESULT_512(2, 3); DECLARE_RESULT_512(3, 3);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x); LOAD_KA_512(1, x); LOAD_KA_512(2, x); LOAD_KA_512(3, x);
					LOAD_KB_512(x, 0); LOAD_KB_512(x, 1); LOAD_KB_512(x, 2); LOAD_KB_512(x, 3);

					MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
					MATMUL_512(0, 2); MATMUL_512(1, 2); MATMUL_512(2, 2); MATMUL_512(3, 2);
					MATMUL_512(0, 3); MATMUL_512(1, 3); MATMUL_512(2, 3); MATMUL_512(3, 3);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x); MASK_LOAD_KA_512(1, x); MASK_LOAD_KA_512(2, x); MASK_LOAD_KA_512(3, x);
					MASK_LOAD_KB_512(x, 0); MASK_LOAD_KB_512(x, 1); MASK_LOAD_KB_512(x, 2); MASK_LOAD_KB_512(x, 3);

					MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
					MATMUL_512(0, 2); MATMUL_512(1, 2); MATMUL_512(2, 2); MATMUL_512(3, 2);
					MATMUL_512(0, 3); MATMUL_512(1, 3); MATMUL_512(2, 3); MATMUL_512(3, 3);
				}
				STORE_REDUCE(0, 0); STORE_REDUCE(1, 0); STORE_REDUCE(2, 0); STORE_REDUCE(3, 0);
				STORE_REDUCE(0, 1); STORE_REDUCE(1, 1); STORE_REDUCE(2, 1); STORE_REDUCE(3, 1);
				STORE_REDUCE(0, 2); STORE_REDUCE(1, 2); STORE_REDUCE(2, 2); STORE_REDUCE(3, 2);
				STORE_REDUCE(0, 3); STORE_REDUCE(1, 3); STORE_REDUCE(2, 3); STORE_REDUCE(3, 3);
			}
			for (; j < n2; j += 2) {
				DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
				DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1); DECLARE_RESULT_512(2, 1); DECLARE_RESULT_512(3, 1);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x); LOAD_KA_512(1, x); LOAD_KA_512(2, x); LOAD_KA_512(3, x);
					LOAD_KB_512(x, 0); LOAD_KB_512(x, 1);

					MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x); MASK_LOAD_KA_512(1, x); MASK_LOAD_KA_512(2, x); MASK_LOAD_KA_512(3, x);
					MASK_LOAD_KB_512(x, 0); MASK_LOAD_KB_512(x, 1);

					MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
				}
				STORE_REDUCE(0, 0); STORE_REDUCE(1, 0); STORE_REDUCE(2, 0); STORE_REDUCE(3, 0);
				STORE_REDUCE(0, 1); STORE_REDUCE(1, 1); STORE_REDUCE(2, 1); STORE_REDUCE(3, 1);

			}
			for (; j < N; j += 1) {
				DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x); LOAD_KA_512(1, x); LOAD_KA_512(2, x); LOAD_KA_512(3, x);
					LOAD_KB_512(x, 0);

					MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x); MASK_LOAD_KA_512(1, x); MASK_LOAD_KA_512(2, x); MASK_LOAD_KA_512(3, x);
					MASK_LOAD_KB_512(x, 0);

					MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
				}
				STORE_REDUCE(0, 0); STORE_REDUCE(1, 0); STORE_REDUCE(2, 0); STORE_REDUCE(3, 0);
			}

		}
		for (; i < m2; i += 2, mi += 2) {
			for (j = 0; j < n4; j += 4) {
				DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0);
				DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1);
				DECLARE_RESULT_512(0, 2); DECLARE_RESULT_512(1, 2);
				DECLARE_RESULT_512(0, 3); DECLARE_RESULT_512(1, 3);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x); LOAD_KA_512(1, x);
					LOAD_KB_512(x, 0); LOAD_KB_512(x, 1); LOAD_KB_512(x, 2); LOAD_KB_512(x, 3);

					MATMUL_512(0, 0); MATMUL_512(1, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1);
					MATMUL_512(0, 2); MATMUL_512(1, 2);
					MATMUL_512(0, 3); MATMUL_512(1, 3);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x); MASK_LOAD_KA_512(1, x);
					MASK_LOAD_KB_512(x, 0); MASK_LOAD_KB_512(x, 1); MASK_LOAD_KB_512(x, 2); MASK_LOAD_KB_512(x, 3);

					MATMUL_512(0, 0); MATMUL_512(1, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1);
					MATMUL_512(0, 2); MATMUL_512(1, 2);
					MATMUL_512(0, 3); MATMUL_512(1, 3);
				}
				STORE_REDUCE(0, 0); STORE_REDUCE(1, 0);
				STORE_REDUCE(0, 1); STORE_REDUCE(1, 1);
				STORE_REDUCE(0, 2); STORE_REDUCE(1, 2);
				STORE_REDUCE(0, 3); STORE_REDUCE(1, 3);
			}
			for (; j < n2; j += 2) {
				DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0);
				DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x); LOAD_KA_512(1, x);
					LOAD_KB_512(x, 0); LOAD_KB_512(x, 1);

					MATMUL_512(0, 0); MATMUL_512(1, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x); MASK_LOAD_KA_512(1, x);
					MASK_LOAD_KB_512(x, 0); MASK_LOAD_KB_512(x, 1);

					MATMUL_512(0, 0); MATMUL_512(1, 0);
					MATMUL_512(0, 1); MATMUL_512(1, 1);
				}
				STORE_REDUCE(0, 0); STORE_REDUCE(1, 0);
				STORE_REDUCE(0, 1); STORE_REDUCE(1, 1);

			}
			for (; j < N; j += 1) {
				DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x); LOAD_KA_512(1, x);
					LOAD_KB_512(x, 0);

					MATMUL_512(0, 0); MATMUL_512(1, 0);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x); MASK_LOAD_KA_512(1, x);
					MASK_LOAD_KB_512(x, 0);

					MATMUL_512(0, 0); MATMUL_512(1, 0);
				}
				STORE_REDUCE(0, 0); STORE_REDUCE(1, 0);
			}
		}
		for (; i < M; i += 1, mi += 1) {
			for (j = 0; j < n4; j += 4) {
				DECLARE_RESULT_512(0, 0);
				DECLARE_RESULT_512(0, 1);
				DECLARE_RESULT_512(0, 2);
				DECLARE_RESULT_512(0, 3);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x);
					LOAD_KB_512(x, 0); LOAD_KB_512(x, 1); LOAD_KB_512(x, 2); LOAD_KB_512(x, 3);

					MATMUL_512(0, 0);
					MATMUL_512(0, 1);
					MATMUL_512(0, 2);
					MATMUL_512(0, 3);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x);
					MASK_LOAD_KB_512(x, 0); MASK_LOAD_KB_512(x, 1); MASK_LOAD_KB_512(x, 2); MASK_LOAD_KB_512(x, 3);


					MATMUL_512(0, 0);
					MATMUL_512(0, 1);
					MATMUL_512(0, 2);
					MATMUL_512(0, 3);
				}
				STORE_REDUCE(0, 0);
				STORE_REDUCE(0, 1);
				STORE_REDUCE(0, 2);
				STORE_REDUCE(0, 3);
			}
			for (; j < n2; j += 2) {
				DECLARE_RESULT_512(0, 0);
				DECLARE_RESULT_512(0, 1);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x);
					LOAD_KB_512(x, 0); LOAD_KB_512(x, 1);

					MATMUL_512(0, 0);
					MATMUL_512(0, 1);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x);
					MASK_LOAD_KB_512(x, 0); MASK_LOAD_KB_512(x, 1);

					MATMUL_512(0, 0);
					MATMUL_512(0, 1);
				}
				STORE_REDUCE(0, 0);
				STORE_REDUCE(0, 1);

			}
			for (; j < N; j += 1) {
				DECLARE_RESULT_512(0, 0);
				for (k = 0; k < k16; k += 16) {
					LOAD_KA_512(0, x);
					LOAD_KB_512(x, 0);

					MATMUL_512(0, 0);
				}
				int remains = K - k;
				if (remains) {
					mask = (1UL << remains) - 1;
					MASK_LOAD_KA_512(0, x);
					MASK_LOAD_KB_512(x, 0);

					MATMUL_512(0, 0);
				}
				STORE_REDUCE(0, 0);
			}
		}
		free(mbuf);
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
