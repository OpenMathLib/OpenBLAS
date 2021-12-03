//  reproduce segfault in zhemv() from zsymv_L_sse2.S
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <sys/mman.h>

#define CALL_ZHEMV zhemv_

void zhemv_(char *UPLO, int *N, double *alpha, double *A, int *LDA,
		double *X, int *INCX, double *beta, double *Y, int *INCY);

int main () {
    
    // zhemv parameters
    char uplo = 'L';
    int n = 14;
    int lda = 16;
    int incx = 1;
    int incy = 1;
    double *A, *X, *Y;    
    double alpha[] = {1, 0};
    double beta[] = {0, 0};

    // other parameters
    int i, j;
    double *data, *data_end, *no_access;
    double real, imag;
    int size;
    size_t len;
    int A_offset;

    size = sizeof(complex double);
    len = lda * lda * size;

    // allocate memory for data
    // use mmap address hints to set up inaccessible memory section following data
    no_access = mmap(NULL, len, PROT_NONE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
    data = mmap(no_access, len, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
    data_end = data + (lda * lda * 2);
    printf("data start/end: %p/%p. Blocked region starts at %p.\n", data, data_end, no_access);

    // set up pointer offsets into data
    A_offset = (lda + 1) * 2;
    A = data + A_offset * 2;			// A starts in the third column of data matrix
    X = data + A_offset + 2;			// X is the second column of data matrix
    Y = (double *)malloc(n * incy * size);	// Y is stored elsewhere
    printf("Address of data: %p; A: %p; X: %p; Y: %p.\n", data, A, X, Y);


    // hermitian matrix
    srand(lda);
    for (j=0; j<lda; j++) {
        real = (double) rand() / RAND_MAX;
        imag = 0;
        data[(j*lda + j) * 2] = real;
        data[(j*lda + j) * 2 + 1] = imag;
        for (i=j+1; i<lda; i++) {
            real = (double) rand() / RAND_MAX;
            imag = (double) rand() / RAND_MAX;
            data[(j*lda + i) * 2] = real;
            data[(j*lda + i) * 2 + 1] = imag;
            data[(i*lda + j) * 2] = real;
            data[(i*lda + j) * 2 + 1] = -imag;
        }
    }

    for (int i=0; i<incy*n*2; i++) {
        Y[i] = 0;
    }

    CALL_ZHEMV(&uplo, &n, alpha, A, &lda, X, &incx, beta, Y, &incy);

    printf("Finished call to zhemv.\n");

    munmap(no_access, len);
    munmap(data, len);

}



