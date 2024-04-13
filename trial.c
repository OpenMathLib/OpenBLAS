#include <stdio.h>
#include <cblas.h>

int main() {
    int n = 4;
    double x[] = {1, 2, 3, 4};
    double y[] = {5, 6, 7, 8};
    double c = 0.5; // cosine of the angle
    double s = 0.5; // sine of the angle
    double a = 2.0; // scalar for axpy and scal

    // Apply the rotation
    cblas_drot(n, x, 1, y, 1, c, s);

    // Perform AXPY operation: y = a*x + y
    cblas_daxpy(n, a, x, 1, y, 1);

    // Swap vectors x and y
    cblas_dswap(n, x, 1, y, 1);

    // Copy vector x to y
    cblas_dcopy(n, x, 1, y, 1);

    // Scale vector x by a
    cblas_dscal(n, a, x, 1);

    // Compute dot product of x and y
    double dot_product = cblas_ddot(n, x, 1, y, 1);

    // Compute sum of absolute values of elements in x
    double asum = cblas_dasum(n, x, 1);

    // Compute Euclidean norm (L2 norm) of vector x
    double norm = cblas_dnrm2(n, x, 1);

    // Find index of maximum absolute value in x
    int index_max = cblas_idamax(n, x, 1);

    // Find index of minimum absolute value in x
    int index_min = cblas_idamin(n, x, 1);

    // Print the results
    printf("Resulting vectors:\n");
    printf("x: ");
    for (int i = 0; i < n; i++) {
        printf("%f ", x[i]);
    }
    printf("\n");
    printf("y: ");
    for (int i = 0; i < n; i++) {
        printf("%f ", y[i]);
    }
    printf("\n");
    printf("Dot product: %f\n", dot_product);
    printf("Asum: %f\n", asum);
    printf("Norm: %f\n", norm);
    printf("Index of max absolute value in x: %d\n", index_max);
    printf("Index of min absolute value in x: %d\n", index_min);

    return 0;
}
