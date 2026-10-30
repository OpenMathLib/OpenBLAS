/***************************************************************************
Copyright (c) 2026, The OpenBLAS Project
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

#include "common.h"
#include <stdint.h>
#if defined(__AVX__)
#include <immintrin.h>
#endif

/*
 * Double precision transposed copy B := alpha * A^T for the OMATCOPY_RT
 * (row-major, ROWM) and OMATCOPY_CT (column-major) kernels.
 *
 * Both orders do the same memory operation, b[j*ldb + i] = alpha * a[i*lda + j],
 * the column-major one with rows and cols swapped. The matrix is walked in
 * 32x32 tiles of 4x4 blocks, transposed in AVX registers.
 */

#define TILE 32

static void dtrans(BLASLONG m, BLASLONG n, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *b, BLASLONG ldb)
{
	BLASLONG i, j, i0, j0, i1, j1;
	BLASLONG m4 = m & ~3, n4 = n & ~3;
#if defined(__AVX__)
	__m256d va = _mm256_set1_pd(alpha);
#endif

	for (i0 = 0; i0 < m4; i0 += TILE) {
		i1 = (i0 + TILE < m4) ? i0 + TILE : m4;
		for (j0 = 0; j0 < n4; j0 += TILE) {
			j1 = (j0 + TILE < n4) ? j0 + TILE : n4;
			for (i = i0; i < i1; i += 4) {
				FLOAT *ap = a + i * lda;
				FLOAT *bp = b + i;
				for (j = j0; j < j1; j += 4) {
#if defined(__AVX__)
					__m256d r0 = _mm256_mul_pd(_mm256_loadu_pd(ap + j), va);
					__m256d r1 = _mm256_mul_pd(_mm256_loadu_pd(ap + lda + j), va);
					__m256d r2 = _mm256_mul_pd(_mm256_loadu_pd(ap + 2 * lda + j), va);
					__m256d r3 = _mm256_mul_pd(_mm256_loadu_pd(ap + 3 * lda + j), va);
					__m256d t0 = _mm256_unpacklo_pd(r0, r1);
					__m256d t1 = _mm256_unpackhi_pd(r0, r1);
					__m256d t2 = _mm256_unpacklo_pd(r2, r3);
					__m256d t3 = _mm256_unpackhi_pd(r2, r3);
					_mm256_storeu_pd(bp + (j + 0) * ldb, _mm256_permute2f128_pd(t0, t2, 0x20));
					_mm256_storeu_pd(bp + (j + 1) * ldb, _mm256_permute2f128_pd(t1, t3, 0x20));
					_mm256_storeu_pd(bp + (j + 2) * ldb, _mm256_permute2f128_pd(t0, t2, 0x31));
					_mm256_storeu_pd(bp + (j + 3) * ldb, _mm256_permute2f128_pd(t1, t3, 0x31));
#else
					BLASLONG ii, jj;
					for (ii = 0; ii < 4; ii++)
						for (jj = 0; jj < 4; jj++)
							bp[(j + jj) * ldb + ii] = alpha * ap[ii * lda + j + jj];
#endif
				}
			}
		}
	}
	for (i = 0; i < m4; i++)
		for (j = n4; j < n; j++)
			b[j * ldb + i] = alpha * a[i * lda + j];
	for (i = m4; i < m; i++)
		for (j = 0; j < n; j++)
			b[j * ldb + i] = alpha * a[i * lda + j];
}

/*
 * Buffers from malloc are often only 16-byte aligned, and with them the
 * blocked loop above runs slower than the generic kernel at some sizes.
 * When a leading dimension is a multiple of 4, every row then has the same
 * alignment, so do the first few rows and columns here and start the
 * blocked part on 32-byte boundaries.
 */
static void dtrans_aligned(BLASLONG m, BLASLONG n, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *b, BLASLONG ldb)
{
	BLASLONG i, j, pi = 0, pj = 0;

	if ((ldb & 3) == 0) pi = ((32 - ((uintptr_t)b & 31)) & 31) / sizeof(FLOAT);
	if ((lda & 3) == 0) pj = ((32 - ((uintptr_t)a & 31)) & 31) / sizeof(FLOAT);
	if (pi > m) pi = m;
	if (pj > n) pj = n;

	for (i = 0; i < pi; i++)
		for (j = 0; j < n; j++)
			b[j * ldb + i] = alpha * a[i * lda + j];
	for (i = pi; i < m; i++)
		for (j = 0; j < pj; j++)
			b[j * ldb + i] = alpha * a[i * lda + j];

	dtrans(m - pi, n - pj, alpha, a + pi * lda + pj, lda, b + pj * ldb + pi, ldb);
}

int CNAME(BLASLONG rows, BLASLONG cols, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *b, BLASLONG ldb)
{
	if (rows <= 0) return 0;
	if (cols <= 0) return 0;
#ifdef ROWM
	dtrans_aligned(rows, cols, alpha, a, lda, b, ldb);
#else
	if (alpha == 0.0) {
		/* same as the generic CT kernel: write zeros without reading A */
		BLASLONG i, j;
		for (i = 0; i < cols; i++)
			for (j = 0; j < rows; j++)
				b[j * ldb + i] = 0.0;
		return 0;
	}
	dtrans_aligned(cols, rows, alpha, a, lda, b, ldb);
#endif
	return 0;
}
