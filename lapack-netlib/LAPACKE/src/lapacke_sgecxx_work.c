/*****************************************************************************
  Copyright (c) 2014, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
* Contents: Native middle-level C interface to LAPACK function sgecxx
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_sgecxx_work( int matrix_layout,
                      char fact, char usesd, lapack_int m, lapack_int n,
                      lapack_int* desel_rows, lapack_int* sel_desel_cols,
                      lapack_int  kmaxfree, float abstol, float reltol,
                      float* a, lapack_int lda, lapack_int* k,
                      float* maxc2nrmk, float* relmaxc2nrmk, float* fnrmk,
                      lapack_int* ipiv, lapack_int* jpiv, float* tau,
                      float* c, lapack_int ldc, float* qrc, lapack_int ldqrc,
                      float* x, lapack_int ldx, float* work, lapack_int lwork,
                      lapack_int* iwork, lapack_int liwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_sgecxx( &fact, &usesd, &m, &n, desel_rows, sel_desel_cols,
                       &kmaxfree, &abstol, &reltol, a, &lda,
                       k, maxc2nrmk, relmaxc2nrmk, fnrmk,
                       ipiv, jpiv, tau, c, &ldc, qrc, &ldqrc, x, &ldx,
                       work, &lwork, iwork, &liwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_logical fact_p = LAPACKE_lsame( fact, 'p' );
        lapack_logical fact_c = LAPACKE_lsame( fact, 'c' );
        lapack_logical fact_x = LAPACKE_lsame( fact, 'x' );
        lapack_int lda_t = MAX(1,m);
        lapack_int ldc_t = MAX(1,m);
        lapack_int ldqrc_t = MAX(1,m);
        lapack_int ldx_t = MAX(1,m);
        float* a_t = NULL;
        float* c_t = NULL;
        float* qrc_t = NULL;
        float* x_t = NULL;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -12;
            LAPACKE_xerbla( "LAPACKE_sgecxx_work", info );
            return info;
        }
        if( ldc < n ) {
            info = -21;
            LAPACKE_xerbla( "LAPACKE_sgecxx_work", info );
            return info;
        }
        if( ldqrc < MIN(m,n) ) {
            info = -23;
            LAPACKE_xerbla( "LAPACKE_sgecxx_work", info );
            return info;
        }
        if( ldx < n ) {
            info = -25;
            LAPACKE_xerbla( "LAPACKE_sgecxx_work", info );
            return info;
        }
        /* Query optimal working array(s) size if requested */
        if( lwork == -1 || liwork == -1 ) {
        LAPACK_sgecxx( &fact, &usesd, &m, &n, desel_rows, sel_desel_cols,
                       &kmaxfree, &abstol, &reltol, a, &lda_t,
                       k, maxc2nrmk, relmaxc2nrmk, fnrmk,
                       ipiv, jpiv, tau, c, &ldc_t, qrc, &ldqrc_t, x, &ldx_t,
                       work, &lwork, iwork, &liwork, &info );
            return (info < 0) ? (info - 1) : info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (float*)LAPACKE_malloc( sizeof(float) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        if( fact_c || fact_x ) {
            c_t = (float*)LAPACKE_malloc( sizeof(float) * ldc_t * MAX(1,n) );
            if( c_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_1;
            }
        }
        if( fact_x ) {
            qrc_t = (float*)LAPACKE_malloc( sizeof(float) * ldqrc_t * MAX(1,MIN(m,n)) );
            if( qrc_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_2;
            }
            x_t = (float*)LAPACKE_malloc( sizeof(float) * ldx_t * MAX(1,n) );
            if( x_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_3;
            }
        }
        /* Transpose input matrices */
        LAPACKE_sge_trans( matrix_layout, m, n, a, lda, a_t, lda_t );
        if( fact_c || fact_x ) {
            LAPACKE_sge_trans( matrix_layout, m, n, c, ldc, c_t, ldc_t );
        }
        if( fact_x ) {
            LAPACKE_sge_trans( matrix_layout, m, MIN(m,n), qrc, ldqrc, qrc_t, ldqrc_t );
            LAPACKE_sge_trans( matrix_layout, m, n, x, ldx, x_t, ldx_t );
        }
        /* Call LAPACK function and adjust info */
        if( fact_p ) {
        LAPACK_sgecxx( &fact, &usesd, &m, &n, desel_rows, sel_desel_cols,
                       &kmaxfree, &abstol, &reltol, a_t, &lda_t,
                       k, maxc2nrmk, relmaxc2nrmk, fnrmk,
                       ipiv, jpiv, tau, c, &ldc_t, qrc, &ldqrc_t, x, &ldx_t,
                       work, &lwork, iwork, &liwork, &info );
        } else if ( fact_c ) {
        LAPACK_sgecxx( &fact, &usesd, &m, &n, desel_rows, sel_desel_cols,
                       &kmaxfree, &abstol, &reltol, a_t, &lda_t,
                       k, maxc2nrmk, relmaxc2nrmk, fnrmk,
                       ipiv, jpiv, tau, c_t, &ldc_t, qrc, &ldqrc_t, x, &ldx_t,
                       work, &lwork, iwork, &liwork, &info );
        } else if ( fact_x ) {
        LAPACK_sgecxx( &fact, &usesd, &m, &n, desel_rows, sel_desel_cols,
                       &kmaxfree, &abstol, &reltol, a_t, &lda_t,
                       k, maxc2nrmk, relmaxc2nrmk, fnrmk,
                       ipiv, jpiv, tau, c_t, &ldc_t, qrc_t, &ldqrc_t, x_t, &ldx_t,
                       work, &lwork, iwork, &liwork, &info );
        }
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, n, a_t, lda_t, a, lda );
        if( fact_c || fact_x ) {
            LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, n, c_t, ldc_t, c, ldc );
        }
        if( fact_x ) {
            LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, MIN(m,n), qrc_t, ldqrc_t, qrc, ldqrc );
            LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, n, x_t, ldx_t, x, ldx );
        }
        /* Release memory and exit */
        if( fact_x ) {
            LAPACKE_free( x_t );
        }
exit_level_3:
        if( fact_x ) {
            LAPACKE_free( qrc_t );
        }
exit_level_2:
        if( fact_c || fact_x ) {
            LAPACKE_free( c_t );
        }
exit_level_1:
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_sgecxx_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_sgecxx_work", info );
    }
    return info;
}
