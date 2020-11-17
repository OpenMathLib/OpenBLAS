#include <stdlib.h>
#include <stdio.h>


#ifdef _OPENMP
    int omp_get_num_procs();
    int omp_get_max_threads();
    void omp_set_num_threads(int);
#endif


#define bint int

void dgemm_(char *transa, char *transb, bint *m, bint *n, bint *k, double *alpha,
    const double *a, bint *lda, const double *b, bint *ldb,
    double *beta, double *c, bint *ldc);

void test(int id)
{
    bint i = 0, m = 1000, n = 800, k = 600, N = m*k + k*n + m*n;
    double alpha = 1.0, beta = 0.0;
    double *A, *B, *C;

    printf("%d\n", id);

    A = (double*)malloc(N*sizeof(double));
    B = A + m*k;
    C = B + k*n;

    for (i = 0; i < N; ++i)
        A[i] = (double) rand() / RAND_MAX;

    dgemm_("N", "N", &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
    dgemm_("T", "N", &m, &n, &k, &alpha, A, &k, B, &k, &beta, C, &m);
    dgemm_("N", "T", &m, &n, &k, &alpha, A, &m, B, &n, &beta, C, &m);
    dgemm_("T", "T", &m, &n, &k, &alpha, A, &k, B, &n, &beta, C, &m);

    free(A);
}

int main()
{
    int i = 0;

#ifdef _OPENMP
    int num_procs = omp_get_num_procs();
    if (num_procs < omp_get_max_threads() * 2)
        omp_set_num_threads(num_procs < 4 ? 1 : num_procs/2);
#endif

    #pragma omp parallel for
    for (i = 0; i < 100; ++i)
        test(i);

    return 0;
}
