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

#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
    #define S0  1
    #define S1 -1
    #define S2  1
    #define S3  1
    #define VFMACC_RR __riscv_vfmsac
    #define VFMACC_RI __riscv_vfmacc
#endif
#if defined(NR) || defined(NC) || defined(TR) || defined(TC)
    #define S0  1
    #define S1  1
    #define S2  1
    #define S3 -1
    #define VFMACC_RR __riscv_vfmacc
    #define VFMACC_RI __riscv_vfmsac
#endif
#if defined(RN) || defined(RT) || defined(CN) || defined(CT)
    #define S0  1
    #define S1  1
    #define S2 -1
    #define S3  1
    #define VFMACC_RR __riscv_vfmacc
    #define VFMACC_RI __riscv_vfnmsac
#endif
#if defined(RR) || defined(RC) || defined(CR) || defined(CC)
    #define S0  1
    #define S1 -1
    #define S2 -1
    #define S3 -1
    #define VFMACC_RR __riscv_vfmsac
    #define VFMACC_RI __riscv_vfnmacc
#endif

#define RISCV_REPEAT_1(INSN, BLEN, ...) \
do { \
    INSN (0, 0, __VA_ARGS__); \
    if (BLEN == 1) break; \
    INSN (1, 1, __VA_ARGS__); \
    if (BLEN == 2) break; \
    INSN (2, 2, __VA_ARGS__); \
    INSN (3, 3, __VA_ARGS__); \
} while (0)

#define RISCV_REPEAT_2(INSN, BLEN, ...) \
do { \
    if (BLEN <= 4) break; \
    INSN (0, 4, __VA_ARGS__); \
    INSN (1, 5, __VA_ARGS__); \
    INSN (2, 6, __VA_ARGS__); \
    INSN (3, 7, __VA_ARGS__); \
} while (0)

#define RISCV_MUL(M, N, DESTr, DESTi, Ar, Ai, Bi, GVL) \
    DESTr##M = RVV_MUL(Ai, Bi[N], GVL); \
    DESTi##M = RVV_MUL(Ar, Bi[N], GVL);

#define RISCV_VFMACC(M, N, DESTr, DESTi, Ar, Ai, Br, GVL) \
    DESTr##M = VFMACC_RR(DESTr##M, Br[N], Ar, GVL); \
    DESTi##M = VFMACC_RI(DESTi##M, Br[N], Ai, GVL);

#define RISCV_ACC_MUL_CONSTR(M, N, DESTr, DESTi, Ar, Ai, B, GVL) \
    DESTr##N = __riscv_vfmacc(DESTr##N, B, Ar##N, GVL); \
    DESTi##N = __riscv_vfmacc(DESTi##N, B, Ai##N, GVL);

#define RISCV_ACC_MUL_CONSTI(M, N, DESTr, DESTi, Ar, Ai, B, GVL) \
    DESTr##N = __riscv_vfnmsac(DESTr##N, B, Ai##N, GVL); \
    DESTi##N = __riscv_vfmacc(DESTi##N, B, Ar##N, GVL);

#define RISCV_LOAD(M, N, DESTr, DESTi, SRC, OFFSET, LDC, GVL) \
    DESTr##N = RVV_LOAD((SRC) + ((OFFSET) + N*(LDC)) * 2, sizeof(FLOAT)*2, GVL); \
    DESTi##N = RVV_LOAD((SRC) + ((OFFSET) + N*(LDC)) * 2 + 1, sizeof(FLOAT)*2, GVL);

#define RISCV_STORE(M, N, DEST, SRCr, SRCi, OFFSET, LDC, GVL) \
    RVV_STORE((DEST) + ((OFFSET) + N*(LDC)) * 2, sizeof(FLOAT)*2, SRCr##N, GVL); \
    RVV_STORE((DEST) + ((OFFSET) + N*(LDC)) * 2 + 1, sizeof(FLOAT)*2, SRCi##N, GVL);

#define RISCV_LOAD_COLUMN(DESTR, DESTI, SRC, OFFSET, GVL) \
    DESTR = RVV_LOAD((SRC) + (OFFSET), sizeof (FLOAT)*2, GVL); \
    DESTI = RVV_LOAD((SRC) + (OFFSET) + 1, sizeof (FLOAT)*2, GVL);

#define COPY_TMP(M, N, DESTr, DESTi, SRCr, SRCi) \
    DESTr##N = SRCr##M; \
    DESTi##N = SRCi##M;

#define RISCV_ADD(M, N, DESTr, DESTi, SRCr, SRCi, GVL) \
    DESTr##N = __riscv_vfadd(DESTr##N, SRCr##M, GVL); \
    DESTi##N = __riscv_vfadd(DESTi##N, SRCi##M, GVL);

#define COPY_ROW(DESTR, DESTI, SRC, OFFSET, LEN) \
do { \
    DESTR[0] = SRC[OFFSET]; \
    DESTI[0] = SRC[OFFSET + 1]; \
    if (LEN == 1) break; \
    DESTR[1] = SRC[OFFSET + 2]; \
    DESTI[1] = SRC[OFFSET + 3]; \
    if (LEN == 2) break; \
    DESTR[2] = SRC[OFFSET + 4]; \
    DESTI[2] = SRC[OFFSET + 5]; \
    DESTR[3] = SRC[OFFSET + 6]; \
    DESTI[3] = SRC[OFFSET + 7]; \
    if (LEN == 4) break; \
    DESTR[4] = SRC[OFFSET + 8]; \
    DESTI[4] = SRC[OFFSET + 9]; \
    DESTR[5] = SRC[OFFSET + 10]; \
    DESTI[5] = SRC[OFFSET + 11]; \
    DESTR[6] = SRC[OFFSET + 12]; \
    DESTI[6] = SRC[OFFSET + 13]; \
    DESTR[7] = SRC[OFFSET + 14]; \
    DESTI[7] = SRC[OFFSET + 15]; \
} while (0)

/* Perform matrix multiplication between submatrices:
   A(m_size,K) * B(K,n_size) = C(m_size,n_size) */

static inline __attribute__((always_inline))
BLASLONG kernel (BLASLONG M, BLASLONG N, BLASLONG K, FLOAT alphar, FLOAT alphai,
                 FLOAT* A, FLOAT* B, FLOAT* C, BLASLONG ldc,
                 BLASLONG m_top, BLASLONG n_top,
                 BLASLONG m_size, BLASLONG n_size)
{
    BLASLONG ai = m_top*K*2;
    BLASLONG bi = n_top*K*2;

    /* b_r[0..n_size-1] = real(B(0, n_top)..B(0, n_top+n_size-1))
       b_i[0..n_size-1] = imag(B(0, n_top)..B(0, n_top+n_size-1)) */
    FLOAT b_r[N_BLOCKSIZE], b_i[N_BLOCKSIZE];
    COPY_ROW (b_r, b_i, B, bi, n_size);
    bi += n_size * 2;

    /* a_r[0..m_size-1] = real(A(m_top, 0)..A(m_top+m_size-1, 0))
       a_i[0..m_size-1] = imag(A(m_top, 0)..A(m_top+m_size-1, 0)) */
    VECTOR_T a_r, a_i;
    RISCV_LOAD_COLUMN (a_r, a_i, A, ai, m_size);
    ai += m_size * 2;

    /* for I = 0..n_size-1
         acc_rI[0..m_size-1] = real(A(m_top..m_top+msize-1, 0) * B(0, ntop+I))
         acc_iI[0..m_size-1] = imag(A(m_top..m_top+msize-1, 0) * B(0, ntop+I)) */
    VECTOR_T tmp_r0, tmp_i0, tmp_r1, tmp_i1, tmp_r2, tmp_i2, tmp_r3, tmp_i3;
    VECTOR_T acc_r0, acc_i0, acc_r1, acc_i1, acc_r2, acc_i2, acc_r3, acc_i3;
    VECTOR_T acc_r4, acc_i4, acc_r5, acc_i5, acc_r6, acc_i6, acc_r7, acc_i7;
    RISCV_REPEAT_1 (RISCV_MUL, n_size, tmp_r, tmp_i, a_r, a_i, b_i, m_size);
    RISCV_REPEAT_1 (RISCV_VFMACC, n_size, tmp_r, tmp_i, a_r, a_i, b_r, m_size);
    RISCV_REPEAT_1 (COPY_TMP, n_size, acc_r, acc_i, tmp_r, tmp_i);
    RISCV_REPEAT_2 (RISCV_MUL, n_size, tmp_r, tmp_i, a_r, a_i, b_i, m_size);
    RISCV_REPEAT_2 (RISCV_VFMACC, n_size, tmp_r, tmp_i, a_r, a_i, b_r, m_size);
    RISCV_REPEAT_2 (COPY_TMP, n_size, acc_r, acc_i, tmp_r, tmp_i);

    for (BLASLONG k = 1; k < K; k++) {
        /* b_r[0..n_size-1] = real(B(k, n_top)..B(k, n_top+n_size-1))
           b_i[0..n_size-1] = imag(B(k, n_top)..B(k, n_top+n_size-1)) */
        COPY_ROW (b_r, b_i, B, bi, n_size);
        bi += n_size * 2;

        /* a_r[0..m_size-1] = real(A(m_top, k)..A(m_top+m_size-1, k))
           a_i[0..m_size-1] = imag(A(m_top, k)..A(m_top+m_size-1, k)) */
        RISCV_LOAD_COLUMN (a_r, a_i, A, ai, m_size);
        ai += m_size * 2;

        /* for I = 0..n_size-1
             acc_rI[0..m_size-1] += real(A(m_top..m_top+msize-1, k) * B(k, ntop+I))
             acc_iI[0..m_size-1] += imag(A(m_top..m_top+msize-1, k) * B(k, ntop+I)) */
        RISCV_REPEAT_1 (RISCV_MUL, n_size, tmp_r, tmp_i, a_r, a_i, b_i, m_size);
        RISCV_REPEAT_1 (RISCV_VFMACC, n_size, tmp_r, tmp_i, a_r, a_i, b_r, m_size);
        RISCV_REPEAT_1 (RISCV_ADD, n_size, acc_r, acc_i, tmp_r, tmp_i, m_size);
        RISCV_REPEAT_2 (RISCV_MUL, n_size, tmp_r, tmp_i, a_r, a_i, b_i, m_size);
        RISCV_REPEAT_2 (RISCV_VFMACC, n_size, tmp_r, tmp_i, a_r, a_i, b_r, m_size);
        RISCV_REPEAT_2 (RISCV_ADD, n_size, acc_r, acc_i, tmp_r, tmp_i, m_size);
    }

    BLASLONG ci = n_top * ldc + m_top;
    VECTOR_T c_r0, c_i0, c_r1, c_i1, c_r2, c_i2, c_r3, c_i3;
    VECTOR_T c_r4, c_i4, c_r5, c_i5, c_r6, c_i6, c_r7, c_i7;

    /* for I = 0..nsize-1
         c_rI[0..m_size-1] = real(C(m_top..m_top+m_size-1, n_top+I))
         c_iI[0..m_size-1] = imag(C(m_top..m_top+m_size-1, n_top+I))
         c_rI[0..m_size-1] += alpha_r * acc_rI[0..m_size-1]
         c_iI[0..m_size-1] += alpha_i * acc_iI[0..m_size-1]
         real(C(mtop..m_top+m_size-1, n_top+I)) = c_rI[0..m_size-1]
         imag(C(mtop..m_top+m_size-1, n_top+I)) = c_iI[0..m_size-1] */
    RISCV_REPEAT_1 (RISCV_LOAD, n_size, c_r, c_i, C, ci, ldc, m_size);
    RISCV_REPEAT_2 (RISCV_LOAD, n_size, c_r, c_i, C, ci, ldc, m_size);
    RISCV_REPEAT_1 (RISCV_ACC_MUL_CONSTR, n_size, c_r, c_i, acc_r, acc_i, alphar, m_size);
    RISCV_REPEAT_2 (RISCV_ACC_MUL_CONSTR, n_size, c_r, c_i, acc_r, acc_i, alphar, m_size);
    RISCV_REPEAT_1 (RISCV_ACC_MUL_CONSTI, n_size, c_r, c_i, acc_r, acc_i, alphai, m_size);
    RISCV_REPEAT_2 (RISCV_ACC_MUL_CONSTI, n_size, c_r, c_i, acc_r, acc_i, alphai, m_size);
    RISCV_REPEAT_1 (RISCV_STORE, n_size, C, c_r, c_i, ci, ldc, m_size);
    RISCV_REPEAT_2 (RISCV_STORE, n_size, C, c_r, c_i, ci, ldc, m_size);

    return m_top + m_size;
}

/* Perform matrix multiplication between submatrices:
   A(M,K) * B(K,n_size) = C(M, n_size) */

static inline __attribute__((always_inline))
BLASLONG kernel_column (BLASLONG M, BLASLONG N, BLASLONG K,
                        FLOAT alphar, FLOAT alphai,
                        FLOAT* A, FLOAT* B, FLOAT* C, BLASLONG ldc,
                        BLASLONG n_top, BLASLONG n_size)
{
    BLASLONG m_top = 0;

    for (BLASLONG i = 0; i < M / M_BLOCKSIZE; i++)
        m_top = kernel (M, N, K, alphar, alphai, A, B, C, ldc, m_top, n_top, M_BLOCKSIZE, n_size);

    if (M & (M_BLOCKSIZE - 1))
        kernel (M, N, K, alphar, alphai, A, B, C, ldc, m_top, n_top, M - m_top, n_size);

    return n_top + n_size;
}

#define xstr(s) str(s)
#define str(s) #s

/* Perform matrix multiplication between matrices:
   A(M,K) * B(K,N) = C(M,N) */

int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT alphar, FLOAT alphai, FLOAT* A, FLOAT* B, FLOAT* C, BLASLONG ldc)
{
    //fprintf(stderr, "%s (with VLV): M=%ld, N=%ld, K=%ld, ldc=%ld\n", xstr(CNAME), M, N, K, ldc);
    BLASLONG n_top = 0;

    for (BLASLONG j = 0; j < N / N_BLOCKSIZE; j++)
        n_top = kernel_column (M, N, K, alphar, alphai, A, B, C, ldc, n_top, N_BLOCKSIZE);

#if N_BLOCKSIZE > 4
    if (N & 4)
        n_top = kernel_column (M, N, K, alphar, alphai, A, B, C, ldc, n_top, 4);
#endif
#if N_BLOCKSIZE > 2
    if (N & 2)
        n_top = kernel_column (M, N, K, alphar, alphai, A, B, C, ldc, n_top, 2);
#endif
#if N_BLOCKSIZE > 1
    if (N & 1)
        kernel_column (M, N, K, alphar, alphai, A, B, C, ldc, n_top, 1);
#endif
    return 0;
}
