#include <cblas.h>

int main(void)
{
    double a[4] = {1, 2, 3, 4}, b[4] = {5, 6, 7, 8}, c[4] = {0};
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0, a, 2, b, 2, 0.0, c, 2);
    return !(c[0] == 19 && c[1] == 22 && c[2] == 43 && c[3] == 50);
}
