#include <string.h>
#include <math.h>
#include <string>

#ifdef __cplusplus
extern "C" {
#endif

void dpotrf_(const char *uplo, const int *n, double *a, const int *lda,
             int *info);

#ifdef __cplusplus
}
#endif

int main(void)
{
	const size_t n = 32768;
	
	// Create positive definite N X N matrix
	double *A = static_cast<double *>(malloc (n * n *sizeof(double)));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (j == i) {
				A[i * n + j] = n - 1;
			} else {
				A[i * n + j] = 1;
			}
		}
	}
	
	// Factorize
	const char *uplo = "U";
	const int int_n = static_cast<int>(n);
	int info;
	dpotrf_(uplo, &int_n, A, &int_n, &info);

	// Free matrix
	free(A);
	
	return 0;
}
