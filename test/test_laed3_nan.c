/*
 * Test for NaN in DSYEVD eigenvectors due to numerical instability
 * in DLAED3 (divide-and-conquer merge) when eigenvalues are nearly degenerate.
 *
 * Build: gcc -o test_laed3_nan test_laed3_nan.c -lopenblas -lm
 * Run:   ./test_laed3_nan
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern void dsyevd_(char *jobz, char *uplo, int *n, double *a, int *lda,
                    double *w, double *work, int *lwork, int *iwork, int *liwork,
                    int *info);

int main(void) {
    int n = 8, lda = 8, lwork = 1 + 6*8 + 2*64, liwork = 3 + 5*8, info;
    double *a = malloc(n*lda*sizeof(double));
    double *w = malloc(n*sizeof(double));
    double *work = malloc(lwork*sizeof(double));
    int *iwork = malloc(liwork*sizeof(int));
    double eps = 1e-8, alpha = 10.0;
    double v[8] = {1.0, 1.0+eps, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i+j*lda] = (i==j ? 1.0 : 0.0) + alpha*v[i]*v[j];
    char jobz = 'V', uplo = 'U';
    dsyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
    if (info) { printf("DSYEVD info=%d\n",info); return 1; }
    int nan=0;
    for (int i=0; i<n; i++) if (isnan(w[i])) nan++;
    for (int i=0; i<n*n; i++) if (isnan(a[i])) nan++;
    if (nan) { printf("FAIL: %d NaN(s)\n", nan); return 1; }
    printf("PASS: No NaN\n");
    free(a); free(w); free(work); free(iwork);
    return 0;
}
