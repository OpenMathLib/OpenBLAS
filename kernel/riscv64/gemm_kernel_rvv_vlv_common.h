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

#define RISCV_REPEAT(INSN, BLEN, ...) \
do { \
    INSN (0, __VA_ARGS__); \
    if (BLEN == 1) break; \
    INSN (1, __VA_ARGS__); \
    if (BLEN == 2) break; \
    INSN (2, __VA_ARGS__); \
    INSN (3, __VA_ARGS__); \
    if (BLEN == 4) break; \
    INSN (4, __VA_ARGS__); \
    INSN (5, __VA_ARGS__); \
    INSN (6, __VA_ARGS__); \
    INSN (7, __VA_ARGS__); \
} while (0)

#define RISCV_MUL(N, DEST, A, B, LEN) \
    DEST##N = RVV_MUL(A, B[N], LEN);

#define RISCV_ACC_MUL(N, DEST, A, B, LEN) \
    DEST##N = RVV_MACC(DEST##N, B[N], A, LEN);

#define RISCV_ACC_MUL_CONST(N, DEST, A, B, LEN) \
    DEST##N = RVV_MACC(DEST##N, B, A##N, LEN);

#define RISCV_LOAD(N, DEST, SRC, OFFSET, LDC, LEN) \
    DEST##N = RVV_LOAD((SRC) + (OFFSET) + N*(LDC), LEN);

#define RISCV_STORE(N, DEST, SRC, OFFSET, LDC, LEN) \
    RVV_STORE((DEST) + (OFFSET) + N*(LDC), SRC##N, LEN);

#define RISCV_LOAD_COLUMN(DEST, SRC, OFFSET, LEN) \
    DEST = RVV_LOAD((SRC) + (OFFSET), LEN);

#define COPY_ROW(DEST, SRC, OFFSET, LEN) \
do { \
    DEST[0] = SRC[OFFSET]; \
    if (LEN == 1) break; \
    DEST[1] = SRC[OFFSET + 1]; \
    if (LEN == 2) break; \
    DEST[2] = SRC[OFFSET + 2]; \
    DEST[3] = SRC[OFFSET + 3]; \
    if (LEN == 4) break; \
    DEST[4] = SRC[OFFSET + 4]; \
    DEST[5] = SRC[OFFSET + 5]; \
    DEST[6] = SRC[OFFSET + 6]; \
    DEST[7] = SRC[OFFSET + 7]; \
} while (0)

/* Perform matrix multiplication between submatrices:
   A(m_size,K) * B(K,n_size) = C(m_size,n_size) */

static inline __attribute__((always_inline))
BLASLONG kernel (BLASLONG M, BLASLONG N, BLASLONG K, FLOAT alpha,
                 FLOAT* A, FLOAT* B, FLOAT* C, BLASLONG ldc,
                 BLASLONG m_top, BLASLONG n_top,
                 BLASLONG m_size, BLASLONG n_size)
{
    BLASLONG ai = m_top*K;
    BLASLONG bi = n_top*K;

    /* b[0..n_size-1] = B(0, n_top)..B(0, n_top+n_size-1) */
    FLOAT b[N_BLOCKSIZE];
    COPY_ROW (b, B, bi, n_size);
    bi += n_size;

    /* a[0..m_size-1] = A(m_top, 0)..A(m_top+m_size-1, 0) */
    VECTOR_T a;
    RISCV_LOAD_COLUMN (a, A, ai, m_size);
    ai += m_size;

    /* for I = 0..n_size-1
         resultI[0..m_size-1] = A(m_top..m_top+msize-1, 0) * B(0, ntop+I) */
    VECTOR_T result0, result1, result2, result3;
    VECTOR_T result4, result5, result6, result7;
    RISCV_REPEAT (RISCV_MUL, n_size, result, a, b, m_size);

    for (BLASLONG k = 1; k < K; k++) {
        /* b[0..n_size-1] = B(k, n_top)..B(k, n_top+n_size-1) */
        COPY_ROW (b, B, bi, n_size);
        bi += n_size;

        /* a[0..m_size-1] = A(m_top, k)..A(m_top+m_size-1, k) */
        RISCV_LOAD_COLUMN (a, A, ai, m_size);
        ai += m_size;

        /* for I = 0..n_size-1
             resultI[0..m_size-1] += A(m_top..m_top+msize-1, k) * B(k, ntop+I) */
        RISCV_REPEAT (RISCV_ACC_MUL, n_size, result, a, b, m_size);
    }

    BLASLONG ci = n_top * ldc + m_top;
    VECTOR_T c0, c1, c2, c3, c4, c5, c6, c7;

    /* for I = 0..nsize-1
         cI[0..m_size-1] = C(m_top..m_top+m_size-1, n_top+I)
         cI[0..m_size-1] += alpha * resultI[0..m_size-1]
         C(mtop..m_top+m_size-1, n_top+I) = cI[0..m_size-1] */
    RISCV_REPEAT (RISCV_LOAD, n_size, c, C, ci, ldc, m_size);
    RISCV_REPEAT (RISCV_ACC_MUL_CONST, n_size, c, result, alpha, m_size);
    RISCV_REPEAT (RISCV_STORE, n_size, C, c, ci, ldc, m_size);

    return m_top + m_size;
}

/* Perform matrix multiplication between submatrices:
   A(M,K) * B(K,n_size) = C(M, n_size) */

static inline __attribute__((always_inline))
BLASLONG kernel_column (BLASLONG M, BLASLONG N, BLASLONG K, FLOAT alpha,
                        FLOAT* A, FLOAT* B, FLOAT* C, BLASLONG ldc,
                        BLASLONG n_top, BLASLONG n_size)
{
    BLASLONG m_top = 0;

    for (BLASLONG i = 0; i < M / M_BLOCKSIZE; i++)
      m_top = kernel (M, N, K, alpha, A, B, C, ldc, m_top, n_top, M_BLOCKSIZE, n_size);

    if (M & (M_BLOCKSIZE - 1))
        kernel (M, N, K, alpha, A, B, C, ldc, m_top, n_top, M - m_top, n_size);

    return n_top + n_size;
}

#define xstr(s) str(s)
#define str(s) #s

/* Perform matrix multiplication between matrices:
   A(M,K) * B(K,N) = C(M,N) */

int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT alpha, FLOAT* A, FLOAT* B, FLOAT* C, BLASLONG ldc)
{
    //fprintf(stderr, "%s (with VLV): M=%ld, N=%ld, K=%ld, ldc=%ld, m_blocksize=%d, n_blocksize=%d\n", xstr(CNAME), M, N, K, ldc, M_BLOCKSIZE, N_BLOCKSIZE);
    BLASLONG n_top = 0;

    for (BLASLONG j = 0; j < N / N_BLOCKSIZE; j++)
        n_top = kernel_column (M, N, K, alpha, A, B, C, ldc, n_top, N_BLOCKSIZE);

#if N_BLOCKSIZE > 4
    if (N & 4)
        n_top = kernel_column (M, N, K, alpha, A, B, C, ldc, n_top, 4);
#endif
#if N_BLOCKSIZE > 2
    if (N & 2)
        n_top = kernel_column (M, N, K, alpha, A, B, C, ldc, n_top, 2);
#endif
#if N_BLOCKSIZE > 1
    if (N & 1)
        kernel_column (M, N, K, alpha, A, B, C, ldc, n_top, 1);
#endif

    return 0;
}
