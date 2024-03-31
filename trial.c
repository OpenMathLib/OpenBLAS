#include <cblas.h>
#include <stdio.h>

int main() {
    int n = 4;
    double x[] = {1, 2, 3, 4};
    double y[] = {5, 6, 7, 8};
    double c = 0.5; // cosine of the angle
    double s = 0.5; // sine of the angle

    // Apply the rotation
    cblas_drot(n, x, 1, y, 1, c, s);

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

    return 0;
}

