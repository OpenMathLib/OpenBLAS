#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>

//------------------------------------------------------------------------------
void fill_rand( int m, int n, double* A, int ld )
{
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            A[ i + j*ld ] = rand() / double(RAND_MAX);
        }
    }
}

//------------------------------------------------------------------------------
inline double max_nan( double x, double y )
{
    return (isnan(y) || (y) >= (x) ? (y) : (x));
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
    int batch_size = 1000;
    int n = 50;
    if (argc > 1)
        batch_size = atoi( argv[1] );
    if (argc > 2)
        n = atoi( argv[2] );
    printf( "batch_size %d, n %d\n", batch_size, n );

    int ld = n;
    double alpha = 3.1416;
    double beta  = 2.7183;

    for (int loops = 0; loops <20; loops ++) {
        
    printf( "init %d\n", loops );
    std::vector<double*> A_array( batch_size ),
                         B_array( batch_size ),
                         C_array( batch_size ),
                         D_array( batch_size );
    for (int i = 0; i < batch_size; ++i) {
        A_array[ i ] = new double[ ld*n ];
        B_array[ i ] = new double[ ld*n ];
        C_array[ i ] = new double[ ld*n ];
        D_array[ i ] = new double[ ld*n ];
        fill_rand( n, n, A_array[ i ], ld );
        fill_rand( n, n, B_array[ i ], ld );
        fill_rand( n, n, C_array[ i ], ld );
        std::copy( C_array[ i ], C_array[ i ] + ld*n, D_array[ i ] );
    }

    printf( "test\n" );
    for (int i = 0; i < batch_size; ++i) {
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
                     alpha, A_array[ i ], ld, B_array[ i ], ld,
                     beta,  C_array[ i ], ld );
    }

    printf( "test OpenMP\n" );
    #pragma omp parallel for
    for (int i = 0; i < batch_size; ++i) {
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
                     alpha, A_array[ i ], ld, B_array[ i ], ld,
                     beta,  D_array[ i ], ld );
    }

    printf( "compare\n" );
    double max_error = 0;
    for (int i = 0; i < batch_size; ++i) {
        // norm( D - C )
        cblas_daxpy( ld*n, -1.0, C_array[ i ], 1, D_array[ i ], 1 );
        double error = cblas_dnrm2( ld*n, D_array[ i ], 1 );
        max_error = max_nan( error, max_error );
    }
    printf( "max error %.2e\n", max_error );

    printf( "delete\n" );
    for (int i = 0; i < batch_size; ++i) {
        delete [] A_array[ i ];
        delete [] B_array[ i ];
        delete [] C_array[ i ];
    }

    printf( "done %d\n", loops );
    } 
    printf( "all done\n");
    return 0;
}

