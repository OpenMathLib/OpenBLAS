typedef long BLASLONG;
typedef float FLOAT;

#ifdef __cplusplus
extern "C" {
#endif

int rdna2_sgemm_f32_host_tuned(float* C, const float* A, const float* Bt, int M, int N, int K);
int rdna2_sgemm_f32_host(float* C, const float* A, const float* Bt, int M, int N, int K);

int sgemm_gfx1201(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, BLASLONG lda,
         FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc);
int sgemm_nn(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, BLASLONG lda,
         FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc);
int sgemm_kernel(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc);

#ifdef __cplusplus
}
#endif

#define CNAME sgemm_gfx1201

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, BLASLONG lda,
         FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
    /* Fallback implementation */
    for (BLASLONG j = 0; j < n; j++) {
        for (BLASLONG i = 0; i < m; i++) {
            FLOAT temp = 0.0f;
            if (beta == 0.0f) {
                for (BLASLONG p = 0; p < k; p++) {
                    temp += a[i + p * lda] * b[p + j * ldb];
                }
                c[i + j * ldc] = alpha * temp;
            } else if (beta == 1.0f) {
                for (BLASLONG p = 0; p < k; p++) {
                    temp += a[i + p * lda] * b[p + j * ldb];
                }
                c[i + j * ldc] += alpha * temp;
            } else {
                for (BLASLONG p = 0; p < k; p++) {
                    temp += a[i + p * lda] * b[p + j * ldb];
                }
                c[i + j * ldc] = beta * c[i + j * ldc] + alpha * temp;
            }
        }
    }
    return 0;
}

int sgemm_nn(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, BLASLONG lda,
         FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
    return CNAME(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

int sgemm_kernel(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {
    /* Try to use our tuned kernel if possible */
    /* For now fallback */
    return CNAME(m, n, k, alpha, a, k, b, n, 0.0f, c, ldc);
}
