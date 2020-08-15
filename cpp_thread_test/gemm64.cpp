#include <iostream>
#include "common.h"
#include "cblas.h"
int main ( int argc, char* argv[] ) {
    const long n = ((long)1 << 31) - 1;
    std::cout << n <<std::endl;
    float*  A = new float[n];
    float*  B = new float[n];
    float*  C = new float[1];
    for(long i =0; i <n; i++){
        A[i] = 1;
        B[i] = 1;

    }
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 1, n,  1.0, A, n, B, 1, 0.0, C, 1);
    std::cout << *C <<std::endl;
    delete[] A;
    delete[] B;
    delete[] C;
    return 0;
}
