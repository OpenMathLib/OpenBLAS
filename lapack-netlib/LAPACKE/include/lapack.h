#ifndef LAPACK_H
#define LAPACK_H

/*
*  Turn on HAVE_LAPACK_CONFIG_H to redefine C-LAPACK datatypes
*/
#ifdef HAVE_LAPACK_CONFIG_H
#include "lapacke_config.h"
#endif

#include "lapacke_mangling.h"

#include <stdlib.h>

/* Complex types are structures equivalent to the
* Fortran complex types COMPLEX(4) and COMPLEX(8).
*
* One can also redefine the types with his own types
* for example by including in the code definitions like
*
* #define lapack_complex_float std::complex<float>
* #define lapack_complex_double std::complex<double>
*
* or define these types in the command line:
*
* -Dlapack_complex_float="std::complex<float>"
* -Dlapack_complex_double="std::complex<double>"
*/

#ifndef LAPACK_COMPLEX_CUSTOM

/* Complex type (single precision) */
#ifndef lapack_complex_float
#ifndef __cplusplus
#include <complex.h>
#else
#include <complex>
#endif
#define lapack_complex_float    float _Complex
#endif

#ifndef lapack_complex_float_real
#define lapack_complex_float_real(z)       (creal(z))
#endif

#ifndef lapack_complex_float_imag
#define lapack_complex_float_imag(z)       (cimag(z))
#endif

/* Complex type (double precision) */
#ifndef lapack_complex_double
#ifndef __cplusplus
#include <complex.h>
#else
#include <complex>
#endif
#define lapack_complex_double   double _Complex
#endif

#ifndef lapack_complex_double_real
#define lapack_complex_double_real(z)      (creal(z))
#endif

#ifndef lapack_complex_double_imag
#define lapack_complex_double_imag(z)       (cimag(z))
#endif

#endif /* LAPACK_COMPLEX_CUSTOM */


#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------*/
#ifndef lapack_int
#define lapack_int     int
#endif

#ifndef lapack_logical
#define lapack_logical lapack_int
#endif

/* f2c, hence clapack and MacOS Accelerate, returns double instead of float
 * for sdot, slange, clange, etc. */
#if defined(LAPACK_F2C)
    typedef double lapack_float_return;
#else
    typedef float lapack_float_return;
#endif


/* Callback logical functions of one, two, or three arguments are used
*  to select eigenvalues to sort to the top left of the Schur form.
*  The value is selected if function returns TRUE (non-zero). */

typedef lapack_logical (*LAPACK_S_SELECT2) ( const float*, const float* );
typedef lapack_logical (*LAPACK_S_SELECT3)
    ( const float*, const float*, const float* );
typedef lapack_logical (*LAPACK_D_SELECT2) ( const double*, const double* );
typedef lapack_logical (*LAPACK_D_SELECT3)
    ( const double*, const double*, const double* );

typedef lapack_logical (*LAPACK_C_SELECT1) ( const lapack_complex_float* );
typedef lapack_logical (*LAPACK_C_SELECT2)
    ( const lapack_complex_float*, const lapack_complex_float* );
typedef lapack_logical (*LAPACK_Z_SELECT1) ( const lapack_complex_double* );
typedef lapack_logical (*LAPACK_Z_SELECT2)
    ( const lapack_complex_double*, const lapack_complex_double* );

#define LAPACK_lsame LAPACK_GLOBAL(lsame,LSAME)
lapack_logical LAPACK_lsame( char* ca,  char* cb,
                              lapack_int lca, lapack_int lcb );


/*----------------------------------------------------------------------------*/
/* This is in alphabetical order (ignoring leading precision). */

#define LAPACK_cbbcsd LAPACK_GLOBAL(cbbcsd,CBBCSD)
void LAPACK_cbbcsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    float* theta,
    float* phi,
    lapack_complex_float* U1, lapack_int const* ldu1,
    lapack_complex_float* U2, lapack_int const* ldu2,
    lapack_complex_float* V1T, lapack_int const* ldv1t,
    lapack_complex_float* V2T, lapack_int const* ldv2t,
    float* B11D,
    float* B11E,
    float* B12D,
    float* B12E,
    float* B21D,
    float* B21E,
    float* B22D,
    float* B22E,
    float* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_dbbcsd LAPACK_GLOBAL(dbbcsd,DBBCSD)
void LAPACK_dbbcsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    double* theta,
    double* phi,
    double* U1, lapack_int const* ldu1,
    double* U2, lapack_int const* ldu2,
    double* V1T, lapack_int const* ldv1t,
    double* V2T, lapack_int const* ldv2t,
    double* B11D,
    double* B11E,
    double* B12D,
    double* B12E,
    double* b21d,
    double* b21e,
    double* b22d,
    double* b22e,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sbbcsd LAPACK_GLOBAL(sbbcsd,SBBCSD)
void LAPACK_sbbcsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    float* theta,
    float* phi,
    float* U1, lapack_int const* ldu1,
    float* U2, lapack_int const* ldu2,
    float* V1T, lapack_int const* ldv1t,
    float* V2T, lapack_int const* ldv2t,
    float* B11D,
    float* B11E,
    float* B12D,
    float* B12E,
    float* B21D,
    float* B21E,
    float* B22D,
    float* B22E,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zbbcsd LAPACK_GLOBAL(zbbcsd,ZBBCSD)
void LAPACK_zbbcsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    double* theta,
    double* phi,
    lapack_complex_double* U1, lapack_int const* ldu1,
    lapack_complex_double* U2, lapack_int const* ldu2,
    lapack_complex_double* V1T, lapack_int const* ldv1t,
    lapack_complex_double* V2T, lapack_int const* ldv2t,
    double* B11D,
    double* B11E,
    double* B12D,
    double* B12E,
    double* B21D,
    double* B21E,
    double* B22D,
    double* B22E,
    double* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_dbdsdc LAPACK_GLOBAL(dbdsdc,DBDSDC)
void LAPACK_dbdsdc(
    char const* uplo, char const* compq,
    lapack_int const* n,
    double* D,
    double* E,
    double* U, lapack_int const* ldu,
    double* VT, lapack_int const* ldvt,
    double* Q, lapack_int* IQ,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sbdsdc LAPACK_GLOBAL(sbdsdc,SBDSDC)
void LAPACK_sbdsdc(
    char const* uplo, char const* compq,
    lapack_int const* n,
    float* D,
    float* E,
    float* U, lapack_int const* ldu,
    float* VT, lapack_int const* ldvt,
    float* Q, lapack_int* IQ,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cbdsqr LAPACK_GLOBAL(cbdsqr,CBDSQR)
void LAPACK_cbdsqr(
    char const* uplo,
    lapack_int const* n, lapack_int const* ncvt, lapack_int const* nru, lapack_int const* ncc,
    float* D,
    float* E,
    lapack_complex_float* VT, lapack_int const* ldvt,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* C, lapack_int const* ldc,
    float* rwork,
    lapack_int* info );

#define LAPACK_dbdsqr LAPACK_GLOBAL(dbdsqr,DBDSQR)
void LAPACK_dbdsqr(
    char const* uplo,
    lapack_int const* n, lapack_int const* ncvt, lapack_int const* nru, lapack_int const* ncc,
    double* D,
    double* E,
    double* VT, lapack_int const* ldvt,
    double* U, lapack_int const* ldu,
    double* C, lapack_int const* ldc,
    double* work,
    lapack_int* info );

#define LAPACK_sbdsqr LAPACK_GLOBAL(sbdsqr,SBDSQR)
void LAPACK_sbdsqr(
    char const* uplo,
    lapack_int const* n, lapack_int const* ncvt, lapack_int const* nru, lapack_int const* ncc,
    float* D,
    float* E,
    float* VT, lapack_int const* ldvt,
    float* U, lapack_int const* ldu,
    float* C, lapack_int const* ldc,
    float* work,
    lapack_int* info );

#define LAPACK_zbdsqr LAPACK_GLOBAL(zbdsqr,ZBDSQR)
void LAPACK_zbdsqr(
    char const* uplo,
    lapack_int const* n, lapack_int const* ncvt, lapack_int const* nru, lapack_int const* ncc,
    double* D,
    double* E,
    lapack_complex_double* VT, lapack_int const* ldvt,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* C, lapack_int const* ldc,
    double* rwork,
    lapack_int* info );

#define LAPACK_dbdsvdx LAPACK_GLOBAL(dbdsvdx,DBDSVDX)
void LAPACK_dbdsvdx(
    char const* uplo, char const* jobz, char const* range,
    lapack_int const* n,
    double const* D,
    double const* E,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* ns,
    double* S,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sbdsvdx LAPACK_GLOBAL(sbdsvdx,SBDSVDX)
void LAPACK_sbdsvdx(
    char const* uplo, char const* jobz, char const* range,
    lapack_int const* n,
    float const* D,
    float const* E,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* ns,
    float* S,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ddisna LAPACK_GLOBAL(ddisna,DDISNA)
void LAPACK_ddisna(
    char const* job,
    lapack_int const* m, lapack_int const* n,
    double const* D,
    double* SEP,
    lapack_int* info );

#define LAPACK_sdisna LAPACK_GLOBAL(sdisna,SDISNA)
void LAPACK_sdisna(
    char const* job,
    lapack_int const* m, lapack_int const* n,
    float const* D,
    float* SEP,
    lapack_int* info );

#define LAPACK_cgbbrd LAPACK_GLOBAL(cgbbrd,CGBBRD)
void LAPACK_cgbbrd(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* ncc, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_float* AB, lapack_int const* ldab,
    float* D,
    float* E,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* PT, lapack_int const* ldpt,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgbbrd LAPACK_GLOBAL(dgbbrd,DGBBRD)
void LAPACK_dgbbrd(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* ncc, lapack_int const* kl, lapack_int const* ku,
    double* AB, lapack_int const* ldab,
    double* D,
    double* E,
    double* Q, lapack_int const* ldq,
    double* PT, lapack_int const* ldpt,
    double* C, lapack_int const* ldc,
    double* work,
    lapack_int* info );

#define LAPACK_sgbbrd LAPACK_GLOBAL(sgbbrd,SGBBRD)
void LAPACK_sgbbrd(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* ncc, lapack_int const* kl, lapack_int const* ku,
    float* AB, lapack_int const* ldab,
    float* D,
    float* E,
    float* Q, lapack_int const* ldq,
    float* PT, lapack_int const* ldpt,
    float* C, lapack_int const* ldc,
    float* work,
    lapack_int* info );

#define LAPACK_zgbbrd LAPACK_GLOBAL(zgbbrd,ZGBBRD)
void LAPACK_zgbbrd(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* ncc, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_double* AB, lapack_int const* ldab,
    double* D,
    double* E,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* PT, lapack_int const* ldpt,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgbcon LAPACK_GLOBAL(cgbcon,CGBCON)
void LAPACK_cgbcon(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_float const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgbcon LAPACK_GLOBAL(dgbcon,DGBCON)
void LAPACK_dgbcon(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    double const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgbcon LAPACK_GLOBAL(sgbcon,SGBCON)
void LAPACK_sgbcon(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    float const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgbcon LAPACK_GLOBAL(zgbcon,ZGBCON)
void LAPACK_zgbcon(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_double const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgbequ LAPACK_GLOBAL(cgbequ,CGBEQU)
void LAPACK_cgbequ(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_dgbequ LAPACK_GLOBAL(dgbequ,DGBEQU)
void LAPACK_dgbequ(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    double const* AB, lapack_int const* ldab,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_sgbequ LAPACK_GLOBAL(sgbequ,SGBEQU)
void LAPACK_sgbequ(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    float const* AB, lapack_int const* ldab,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_zgbequ LAPACK_GLOBAL(zgbequ,ZGBEQU)
void LAPACK_zgbequ(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_cgbequb LAPACK_GLOBAL(cgbequb,CGBEQUB)
void LAPACK_cgbequb(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_dgbequb LAPACK_GLOBAL(dgbequb,DGBEQUB)
void LAPACK_dgbequb(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    double const* AB, lapack_int const* ldab,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_sgbequb LAPACK_GLOBAL(sgbequb,SGBEQUB)
void LAPACK_sgbequb(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    float const* AB, lapack_int const* ldab,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_zgbequb LAPACK_GLOBAL(zgbequb,ZGBEQUB)
void LAPACK_zgbequb(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_cgbrfs LAPACK_GLOBAL(cgbrfs,CGBRFS)
void LAPACK_cgbrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_float const* AB, lapack_int const* ldab,
    lapack_complex_float const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgbrfs LAPACK_GLOBAL(dgbrfs,DGBRFS)
void LAPACK_dgbrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab,
    double const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgbrfs LAPACK_GLOBAL(sgbrfs,SGBRFS)
void LAPACK_sgbrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    float const* AB, lapack_int const* ldab,
    float const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgbrfs LAPACK_GLOBAL(zgbrfs,ZGBRFS)
void LAPACK_zgbrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_double const* AB, lapack_int const* ldab,
    lapack_complex_double const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgbrfsx LAPACK_GLOBAL(cgbrfsx,CGBRFSX)
void LAPACK_cgbrfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_float const* AB, lapack_int const* ldab,
    lapack_complex_float const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    float* R,
    float* C,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgbrfsx LAPACK_GLOBAL(dgbrfsx,DGBRFSX)
void LAPACK_dgbrfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab,
    double const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    double* R,
    double* C,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgbrfsx LAPACK_GLOBAL(sgbrfsx,SGBRFSX)
void LAPACK_sgbrfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    float const* AB, lapack_int const* ldab,
    float const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    float* R,
    float* C,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgbrfsx LAPACK_GLOBAL(zgbrfsx,ZGBRFSX)
void LAPACK_zgbrfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_double const* AB, lapack_int const* ldab,
    lapack_complex_double const* AFB, lapack_int const* ldafb, lapack_int const* ipiv,
    double* R,
    double* C,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgbsv LAPACK_GLOBAL(cgbsv,CGBSV)
void LAPACK_cgbsv(
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_float* AB, lapack_int const* ldab, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dgbsv LAPACK_GLOBAL(dgbsv,DGBSV)
void LAPACK_dgbsv(
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    double* AB, lapack_int const* ldab, lapack_int* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sgbsv LAPACK_GLOBAL(sgbsv,SGBSV)
void LAPACK_sgbsv(
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    float* AB, lapack_int const* ldab, lapack_int* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zgbsv LAPACK_GLOBAL(zgbsv,ZGBSV)
void LAPACK_zgbsv(
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_double* AB, lapack_int const* ldab, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cgbsvx LAPACK_GLOBAL(cgbsvx,CGBSVX)
void LAPACK_cgbsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgbsvx LAPACK_GLOBAL(dgbsvx,DGBSVX)
void LAPACK_dgbsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    double* AB, lapack_int const* ldab,
    double* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgbsvx LAPACK_GLOBAL(sgbsvx,SGBSVX)
void LAPACK_sgbsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    float* AB, lapack_int const* ldab,
    float* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgbsvx LAPACK_GLOBAL(zgbsvx,ZGBSVX)
void LAPACK_zgbsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgbsvxx LAPACK_GLOBAL(cgbsvxx,CGBSVXX)
void LAPACK_cgbsvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgbsvxx LAPACK_GLOBAL(dgbsvxx,DGBSVXX)
void LAPACK_dgbsvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    double* AB, lapack_int const* ldab,
    double* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgbsvxx LAPACK_GLOBAL(sgbsvxx,SGBSVXX)
void LAPACK_sgbsvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    float* AB, lapack_int const* ldab,
    float* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgbsvxx LAPACK_GLOBAL(zgbsvxx,ZGBSVXX)
void LAPACK_zgbsvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* AFB, lapack_int const* ldafb, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgbtrf LAPACK_GLOBAL(cgbtrf,CGBTRF)
void LAPACK_cgbtrf(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_float* AB, lapack_int const* ldab, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_dgbtrf LAPACK_GLOBAL(dgbtrf,DGBTRF)
void LAPACK_dgbtrf(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    double* AB, lapack_int const* ldab, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_sgbtrf LAPACK_GLOBAL(sgbtrf,SGBTRF)
void LAPACK_sgbtrf(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    float* AB, lapack_int const* ldab, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_zgbtrf LAPACK_GLOBAL(zgbtrf,ZGBTRF)
void LAPACK_zgbtrf(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_double* AB, lapack_int const* ldab, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_cgbtrs LAPACK_GLOBAL(cgbtrs,CGBTRS)
void LAPACK_cgbtrs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_float const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dgbtrs LAPACK_GLOBAL(dgbtrs,DGBTRS)
void LAPACK_dgbtrs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sgbtrs LAPACK_GLOBAL(sgbtrs,SGBTRS)
void LAPACK_sgbtrs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    float const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zgbtrs LAPACK_GLOBAL(zgbtrs,ZGBTRS)
void LAPACK_zgbtrs(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    lapack_complex_double const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cgebak LAPACK_GLOBAL(cgebak,CGEBAK)
void LAPACK_cgebak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float const* scale, lapack_int const* m,
    lapack_complex_float* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_dgebak LAPACK_GLOBAL(dgebak,DGEBAK)
void LAPACK_dgebak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double const* scale, lapack_int const* m,
    double* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_sgebak LAPACK_GLOBAL(sgebak,SGEBAK)
void LAPACK_sgebak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float const* scale, lapack_int const* m,
    float* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_zgebak LAPACK_GLOBAL(zgebak,ZGEBAK)
void LAPACK_zgebak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double const* scale, lapack_int const* m,
    lapack_complex_double* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_cgebal LAPACK_GLOBAL(cgebal,CGEBAL)
void LAPACK_cgebal(
    char const* job,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ilo, lapack_int* ihi,
    float* scale,
    lapack_int* info );

#define LAPACK_dgebal LAPACK_GLOBAL(dgebal,DGEBAL)
void LAPACK_dgebal(
    char const* job,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* ilo, lapack_int* ihi,
    double* scale,
    lapack_int* info );

#define LAPACK_sgebal LAPACK_GLOBAL(sgebal,SGEBAL)
void LAPACK_sgebal(
    char const* job,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* ilo, lapack_int* ihi,
    float* scale,
    lapack_int* info );

#define LAPACK_zgebal LAPACK_GLOBAL(zgebal,ZGEBAL)
void LAPACK_zgebal(
    char const* job,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ilo, lapack_int* ihi,
    double* scale,
    lapack_int* info );

#define LAPACK_cgebrd LAPACK_GLOBAL(cgebrd,CGEBRD)
void LAPACK_cgebrd(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* D,
    float* E,
    lapack_complex_float* tauq,
    lapack_complex_float* taup,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgebrd LAPACK_GLOBAL(dgebrd,DGEBRD)
void LAPACK_dgebrd(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* D,
    double* E,
    double* tauq,
    double* taup,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgebrd LAPACK_GLOBAL(sgebrd,SGEBRD)
void LAPACK_sgebrd(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* D,
    float* E,
    float* tauq,
    float* taup,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgebrd LAPACK_GLOBAL(zgebrd,ZGEBRD)
void LAPACK_zgebrd(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* D,
    double* E,
    lapack_complex_double* tauq,
    lapack_complex_double* taup,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgecon LAPACK_GLOBAL(cgecon,CGECON)
void LAPACK_cgecon(
    char const* norm,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgecon LAPACK_GLOBAL(dgecon,DGECON)
void LAPACK_dgecon(
    char const* norm,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgecon LAPACK_GLOBAL(sgecon,SGECON)
void LAPACK_sgecon(
    char const* norm,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgecon LAPACK_GLOBAL(zgecon,ZGECON)
void LAPACK_zgecon(
    char const* norm,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgeequ LAPACK_GLOBAL(cgeequ,CGEEQU)
void LAPACK_cgeequ(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_dgeequ LAPACK_GLOBAL(dgeequ,DGEEQU)
void LAPACK_dgeequ(
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_sgeequ LAPACK_GLOBAL(sgeequ,SGEEQU)
void LAPACK_sgeequ(
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_zgeequ LAPACK_GLOBAL(zgeequ,ZGEEQU)
void LAPACK_zgeequ(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_cgeequb LAPACK_GLOBAL(cgeequb,CGEEQUB)
void LAPACK_cgeequb(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_dgeequb LAPACK_GLOBAL(dgeequb,DGEEQUB)
void LAPACK_dgeequb(
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_sgeequb LAPACK_GLOBAL(sgeequb,SGEEQUB)
void LAPACK_sgeequb(
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* R,
    float* C,
    float* rowcnd,
    float* colcnd,
    float* amax,
    lapack_int* info );

#define LAPACK_zgeequb LAPACK_GLOBAL(zgeequb,ZGEEQUB)
void LAPACK_zgeequb(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* R,
    double* C,
    double* rowcnd,
    double* colcnd,
    double* amax,
    lapack_int* info );

#define LAPACK_cgees LAPACK_GLOBAL(cgees,CGEES)
void LAPACK_cgees(
    char const* jobvs, char const* sort, LAPACK_C_SELECT1 select,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* sdim,
    lapack_complex_float* W,
    lapack_complex_float* VS, lapack_int const* ldvs,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_dgees LAPACK_GLOBAL(dgees,DGEES)
void LAPACK_dgees(
    char const* jobvs, char const* sort, LAPACK_D_SELECT2 select,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* sdim,
    double* WR,
    double* WI,
    double* VS, lapack_int const* ldvs,
    double* work, lapack_int const* lwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_sgees LAPACK_GLOBAL(sgees,SGEES)
void LAPACK_sgees(
    char const* jobvs, char const* sort, LAPACK_S_SELECT2 select,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* sdim,
    float* WR,
    float* WI,
    float* VS, lapack_int const* ldvs,
    float* work, lapack_int const* lwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_zgees LAPACK_GLOBAL(zgees,ZGEES)
void LAPACK_zgees(
    char const* jobvs, char const* sort, LAPACK_Z_SELECT1 select,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* sdim,
    lapack_complex_double* W,
    lapack_complex_double* VS, lapack_int const* ldvs,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_cgeesx LAPACK_GLOBAL(cgeesx,CGEESX)
void LAPACK_cgeesx(
    char const* jobvs, char const* sort, LAPACK_C_SELECT1 select, char const* sense,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* sdim,
    lapack_complex_float* W,
    lapack_complex_float* VS, lapack_int const* ldvs,
    float* rconde,
    float* rcondv,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_dgeesx LAPACK_GLOBAL(dgeesx,DGEESX)
void LAPACK_dgeesx(
    char const* jobvs, char const* sort, LAPACK_D_SELECT2 select, char const* sense,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* sdim,
    double* WR,
    double* WI,
    double* VS, lapack_int const* ldvs,
    double* rconde,
    double* rcondv,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_sgeesx LAPACK_GLOBAL(sgeesx,SGEESX)
void LAPACK_sgeesx(
    char const* jobvs, char const* sort, LAPACK_S_SELECT2 select, char const* sense,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* sdim,
    float* WR,
    float* WI,
    float* VS, lapack_int const* ldvs,
    float* rconde,
    float* rcondv,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_zgeesx LAPACK_GLOBAL(zgeesx,ZGEESX)
void LAPACK_zgeesx(
    char const* jobvs, char const* sort, LAPACK_Z_SELECT1 select, char const* sense,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* sdim,
    lapack_complex_double* W,
    lapack_complex_double* VS, lapack_int const* ldvs,
    double* rconde,
    double* rcondv,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_cgeev LAPACK_GLOBAL(cgeev,CGEEV)
void LAPACK_cgeev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* W,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgeev LAPACK_GLOBAL(dgeev,DGEEV)
void LAPACK_dgeev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* WR,
    double* WI,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgeev LAPACK_GLOBAL(sgeev,SGEEV)
void LAPACK_sgeev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* WR,
    float* WI,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgeev LAPACK_GLOBAL(zgeev,ZGEEV)
void LAPACK_zgeev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* W,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgeevx LAPACK_GLOBAL(cgeevx,CGEEVX)
void LAPACK_cgeevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* W,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    float* scale,
    float* abnrm,
    float* rconde,
    float* rcondv,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgeevx LAPACK_GLOBAL(dgeevx,DGEEVX)
void LAPACK_dgeevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* WR,
    double* WI,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    double* scale,
    double* abnrm,
    double* rconde,
    double* rcondv,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgeevx LAPACK_GLOBAL(sgeevx,SGEEVX)
void LAPACK_sgeevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* WR,
    float* WI,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    float* scale,
    float* abnrm,
    float* rconde,
    float* rcondv,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgeevx LAPACK_GLOBAL(zgeevx,ZGEEVX)
void LAPACK_zgeevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* W,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    double* scale,
    double* abnrm,
    double* rconde,
    double* rcondv,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgehrd LAPACK_GLOBAL(cgehrd,CGEHRD)
void LAPACK_cgehrd(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgehrd LAPACK_GLOBAL(dgehrd,DGEHRD)
void LAPACK_dgehrd(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double* A, lapack_int const* lda,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgehrd LAPACK_GLOBAL(sgehrd,SGEHRD)
void LAPACK_sgehrd(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float* A, lapack_int const* lda,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgehrd LAPACK_GLOBAL(zgehrd,ZGEHRD)
void LAPACK_zgehrd(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgejsv LAPACK_GLOBAL(cgejsv,CGEJSV)
void LAPACK_cgejsv(
    char const* joba, char const* jobu, char const* jobv, char const* jobr, char const* jobt, char const* jobp,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* SVA,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* V, lapack_int const* ldv,
    lapack_complex_float* cwork, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dgejsv LAPACK_GLOBAL(dgejsv,DGEJSV)
void LAPACK_dgejsv(
    char const* joba, char const* jobu, char const* jobv, char const* jobr, char const* jobt, char const* jobp,
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* SVA,
    double* U, lapack_int const* ldu,
    double* V, lapack_int const* ldv,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgejsv LAPACK_GLOBAL(sgejsv,SGEJSV)
void LAPACK_sgejsv(
    char const* joba, char const* jobu, char const* jobv, char const* jobr, char const* jobt, char const* jobp,
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* SVA,
    float* U, lapack_int const* ldu,
    float* V, lapack_int const* ldv,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgejsv LAPACK_GLOBAL(zgejsv,ZGEJSV)
void LAPACK_zgejsv(
    char const* joba, char const* jobu, char const* jobv, char const* jobr, char const* jobt, char const* jobp,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* SVA,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* V, lapack_int const* ldv,
    lapack_complex_double* cwork, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cgelq LAPACK_GLOBAL(cgelq,CGELQ)
void LAPACK_cgelq(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* T, lapack_int const* tsize,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgelq LAPACK_GLOBAL(dgelq,DGELQ)
void LAPACK_dgelq(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* T, lapack_int const* tsize,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgelq LAPACK_GLOBAL(sgelq,SGELQ)
void LAPACK_sgelq(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* T, lapack_int const* tsize,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgelq LAPACK_GLOBAL(zgelq,ZGELQ)
void LAPACK_zgelq(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* T, lapack_int const* tsize,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgelq2 LAPACK_GLOBAL(cgelq2,CGELQ2)
void LAPACK_cgelq2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dgelq2 LAPACK_GLOBAL(dgelq2,DGELQ2)
void LAPACK_dgelq2(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work,
    lapack_int* info );

#define LAPACK_sgelq2 LAPACK_GLOBAL(sgelq2,SGELQ2)
void LAPACK_sgelq2(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work,
    lapack_int* info );

#define LAPACK_zgelq2 LAPACK_GLOBAL(zgelq2,ZGELQ2)
void LAPACK_zgelq2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cgelqf LAPACK_GLOBAL(cgelqf,CGELQF)
void LAPACK_cgelqf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgelqf LAPACK_GLOBAL(dgelqf,DGELQF)
void LAPACK_dgelqf(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgelqf LAPACK_GLOBAL(sgelqf,SGELQF)
void LAPACK_sgelqf(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgelqf LAPACK_GLOBAL(zgelqf,ZGELQF)
void LAPACK_zgelqf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgels LAPACK_GLOBAL(cgels,CGELS)
void LAPACK_cgels(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgels LAPACK_GLOBAL(dgels,DGELS)
void LAPACK_dgels(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgels LAPACK_GLOBAL(sgels,SGELS)
void LAPACK_sgels(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgels LAPACK_GLOBAL(zgels,ZGELS)
void LAPACK_zgels(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgelsd LAPACK_GLOBAL(cgelsd,CGELSD)
void LAPACK_cgelsd(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float* S,
    float const* rcond, lapack_int* rank,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dgelsd LAPACK_GLOBAL(dgelsd,DGELSD)
void LAPACK_dgelsd(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* S,
    double const* rcond, lapack_int* rank,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgelsd LAPACK_GLOBAL(sgelsd,SGELSD)
void LAPACK_sgelsd(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* S,
    float const* rcond, lapack_int* rank,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgelsd LAPACK_GLOBAL(zgelsd,ZGELSD)
void LAPACK_zgelsd(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double* S,
    double const* rcond, lapack_int* rank,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cgelss LAPACK_GLOBAL(cgelss,CGELSS)
void LAPACK_cgelss(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float* S,
    float const* rcond, lapack_int* rank,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgelss LAPACK_GLOBAL(dgelss,DGELSS)
void LAPACK_dgelss(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* S,
    double const* rcond, lapack_int* rank,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgelss LAPACK_GLOBAL(sgelss,SGELSS)
void LAPACK_sgelss(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* S,
    float const* rcond, lapack_int* rank,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgelss LAPACK_GLOBAL(zgelss,ZGELSS)
void LAPACK_zgelss(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double* S,
    double const* rcond, lapack_int* rank,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgelsy LAPACK_GLOBAL(cgelsy,CGELSY)
void LAPACK_cgelsy(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb, lapack_int* JPVT,
    float const* rcond, lapack_int* rank,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgelsy LAPACK_GLOBAL(dgelsy,DGELSY)
void LAPACK_dgelsy(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb, lapack_int* JPVT,
    double const* rcond, lapack_int* rank,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgelsy LAPACK_GLOBAL(sgelsy,SGELSY)
void LAPACK_sgelsy(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb, lapack_int* JPVT,
    float const* rcond, lapack_int* rank,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgelsy LAPACK_GLOBAL(zgelsy,ZGELSY)
void LAPACK_zgelsy(
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb, lapack_int* JPVT,
    double const* rcond, lapack_int* rank,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgemlq LAPACK_GLOBAL(cgemlq,CGEMLQ)
void LAPACK_cgemlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* T, lapack_int const* tsize,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgemlq LAPACK_GLOBAL(dgemlq,DGEMLQ)
void LAPACK_dgemlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* A, lapack_int const* lda,
    double const* T, lapack_int const* tsize,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgemlq LAPACK_GLOBAL(sgemlq,SGEMLQ)
void LAPACK_sgemlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* A, lapack_int const* lda,
    float const* T, lapack_int const* tsize,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgemlq LAPACK_GLOBAL(zgemlq,ZGEMLQ)
void LAPACK_zgemlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* T, lapack_int const* tsize,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgemqr LAPACK_GLOBAL(cgemqr,CGEMQR)
void LAPACK_cgemqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* T, lapack_int const* tsize,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgemqr LAPACK_GLOBAL(dgemqr,DGEMQR)
void LAPACK_dgemqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* A, lapack_int const* lda,
    double const* T, lapack_int const* tsize,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgemqr LAPACK_GLOBAL(sgemqr,SGEMQR)
void LAPACK_sgemqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* A, lapack_int const* lda,
    float const* T, lapack_int const* tsize,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgemqr LAPACK_GLOBAL(zgemqr,ZGEMQR)
void LAPACK_zgemqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* T, lapack_int const* tsize,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgemqrt LAPACK_GLOBAL(cgemqrt,CGEMQRT)
void LAPACK_cgemqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* nb,
    lapack_complex_float const* V, lapack_int const* ldv,
    lapack_complex_float const* T, lapack_int const* ldt,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dgemqrt LAPACK_GLOBAL(dgemqrt,DGEMQRT)
void LAPACK_dgemqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* nb,
    double const* V, lapack_int const* ldv,
    double const* T, lapack_int const* ldt,
    double* C, lapack_int const* ldc,
    double* work,
    lapack_int* info );

#define LAPACK_sgemqrt LAPACK_GLOBAL(sgemqrt,SGEMQRT)
void LAPACK_sgemqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* nb,
    float const* V, lapack_int const* ldv,
    float const* T, lapack_int const* ldt,
    float* C, lapack_int const* ldc,
    float* work,
    lapack_int* info );

#define LAPACK_zgemqrt LAPACK_GLOBAL(zgemqrt,ZGEMQRT)
void LAPACK_zgemqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* nb,
    lapack_complex_double const* V, lapack_int const* ldv,
    lapack_complex_double const* T, lapack_int const* ldt,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cgeql2 LAPACK_GLOBAL(cgeql2,CGEQL2)
void LAPACK_cgeql2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dgeql2 LAPACK_GLOBAL(dgeql2,DGEQL2)
void LAPACK_dgeql2(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work,
    lapack_int* info );

#define LAPACK_sgeql2 LAPACK_GLOBAL(sgeql2,SGEQL2)
void LAPACK_sgeql2(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work,
    lapack_int* info );

#define LAPACK_zgeql2 LAPACK_GLOBAL(zgeql2,ZGEQL2)
void LAPACK_zgeql2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cgeqlf LAPACK_GLOBAL(cgeqlf,CGEQLF)
void LAPACK_cgeqlf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgeqlf LAPACK_GLOBAL(dgeqlf,DGEQLF)
void LAPACK_dgeqlf(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgeqlf LAPACK_GLOBAL(sgeqlf,SGEQLF)
void LAPACK_sgeqlf(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgeqlf LAPACK_GLOBAL(zgeqlf,ZGEQLF)
void LAPACK_zgeqlf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgeqpf LAPACK_GLOBAL(sgeqpf,SGEQPF)
void LAPACK_sgeqpf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
                    lapack_int* jpvt, float* tau, float* work,
                    lapack_int *info );

#define LAPACK_dgeqpf LAPACK_GLOBAL(dgeqpf,DGEQPF)
void LAPACK_dgeqpf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
                    lapack_int* jpvt, double* tau, double* work,
                    lapack_int *info );

#define LAPACK_cgeqpf LAPACK_GLOBAL(cgeqpf,CGEQPF)
void LAPACK_cgeqpf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
                    lapack_int* lda, lapack_int* jpvt,
                    lapack_complex_float* tau, lapack_complex_float* work,
                    float* rwork, lapack_int *info );

#define LAPACK_zgeqpf LAPACK_GLOBAL(zgeqpf,ZGEQPF)
void LAPACK_zgeqpf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
                    lapack_int* lda, lapack_int* jpvt,
                    lapack_complex_double* tau, lapack_complex_double* work,
                    double* rwork, lapack_int *info );

#define LAPACK_cgeqp3 LAPACK_GLOBAL(cgeqp3,CGEQP3)
void LAPACK_cgeqp3(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* JPVT,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgeqp3 LAPACK_GLOBAL(dgeqp3,DGEQP3)
void LAPACK_dgeqp3(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* JPVT,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgeqp3 LAPACK_GLOBAL(sgeqp3,SGEQP3)
void LAPACK_sgeqp3(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* JPVT,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgeqp3 LAPACK_GLOBAL(zgeqp3,ZGEQP3)
void LAPACK_zgeqp3(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* JPVT,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgeqr LAPACK_GLOBAL(cgeqr,CGEQR)
void LAPACK_cgeqr(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* T, lapack_int const* tsize,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgeqr LAPACK_GLOBAL(dgeqr,DGEQR)
void LAPACK_dgeqr(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* T, lapack_int const* tsize,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgeqr LAPACK_GLOBAL(sgeqr,SGEQR)
void LAPACK_sgeqr(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* T, lapack_int const* tsize,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgeqr LAPACK_GLOBAL(zgeqr,ZGEQR)
void LAPACK_zgeqr(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* T, lapack_int const* tsize,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgeqr2 LAPACK_GLOBAL(cgeqr2,CGEQR2)
void LAPACK_cgeqr2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dgeqr2 LAPACK_GLOBAL(dgeqr2,DGEQR2)
void LAPACK_dgeqr2(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work,
    lapack_int* info );

#define LAPACK_sgeqr2 LAPACK_GLOBAL(sgeqr2,SGEQR2)
void LAPACK_sgeqr2(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work,
    lapack_int* info );

#define LAPACK_zgeqr2 LAPACK_GLOBAL(zgeqr2,ZGEQR2)
void LAPACK_zgeqr2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cgeqrf LAPACK_GLOBAL(cgeqrf,CGEQRF)
void LAPACK_cgeqrf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgeqrf LAPACK_GLOBAL(dgeqrf,DGEQRF)
void LAPACK_dgeqrf(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgeqrf LAPACK_GLOBAL(sgeqrf,SGEQRF)
void LAPACK_sgeqrf(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgeqrf LAPACK_GLOBAL(zgeqrf,ZGEQRF)
void LAPACK_zgeqrf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgeqrfp LAPACK_GLOBAL(cgeqrfp,CGEQRFP)
void LAPACK_cgeqrfp(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgeqrfp LAPACK_GLOBAL(dgeqrfp,DGEQRFP)
void LAPACK_dgeqrfp(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgeqrfp LAPACK_GLOBAL(sgeqrfp,SGEQRFP)
void LAPACK_sgeqrfp(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgeqrfp LAPACK_GLOBAL(zgeqrfp,ZGEQRFP)
void LAPACK_zgeqrfp(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgeqrt LAPACK_GLOBAL(cgeqrt,CGEQRT)
void LAPACK_cgeqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* nb,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dgeqrt LAPACK_GLOBAL(dgeqrt,DGEQRT)
void LAPACK_dgeqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* nb,
    double* A, lapack_int const* lda,
    double* T, lapack_int const* ldt,
    double* work,
    lapack_int* info );

#define LAPACK_sgeqrt LAPACK_GLOBAL(sgeqrt,SGEQRT)
void LAPACK_sgeqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* nb,
    float* A, lapack_int const* lda,
    float* T, lapack_int const* ldt,
    float* work,
    lapack_int* info );

#define LAPACK_zgeqrt LAPACK_GLOBAL(zgeqrt,ZGEQRT)
void LAPACK_zgeqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* nb,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cgeqrt2 LAPACK_GLOBAL(cgeqrt2,CGEQRT2)
void LAPACK_cgeqrt2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_dgeqrt2 LAPACK_GLOBAL(dgeqrt2,DGEQRT2)
void LAPACK_dgeqrt2(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_sgeqrt2 LAPACK_GLOBAL(sgeqrt2,SGEQRT2)
void LAPACK_sgeqrt2(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_zgeqrt2 LAPACK_GLOBAL(zgeqrt2,ZGEQRT2)
void LAPACK_zgeqrt2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_cgeqrt3 LAPACK_GLOBAL(cgeqrt3,CGEQRT3)
void LAPACK_cgeqrt3(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_dgeqrt3 LAPACK_GLOBAL(dgeqrt3,DGEQRT3)
void LAPACK_dgeqrt3(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_sgeqrt3 LAPACK_GLOBAL(sgeqrt3,SGEQRT3)
void LAPACK_sgeqrt3(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_zgeqrt3 LAPACK_GLOBAL(zgeqrt3,ZGEQRT3)
void LAPACK_zgeqrt3(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_cgerfs LAPACK_GLOBAL(cgerfs,CGERFS)
void LAPACK_cgerfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgerfs LAPACK_GLOBAL(dgerfs,DGERFS)
void LAPACK_dgerfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgerfs LAPACK_GLOBAL(sgerfs,SGERFS)
void LAPACK_sgerfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgerfs LAPACK_GLOBAL(zgerfs,ZGERFS)
void LAPACK_zgerfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgerfsx LAPACK_GLOBAL(cgerfsx,CGERFSX)
void LAPACK_cgerfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    float const* R,
    float const* C,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgerfsx LAPACK_GLOBAL(dgerfsx,DGERFSX)
void LAPACK_dgerfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    double const* R,
    double const* C,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgerfsx LAPACK_GLOBAL(sgerfsx,SGERFSX)
void LAPACK_sgerfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    float const* R,
    float const* C,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgerfsx LAPACK_GLOBAL(zgerfsx,ZGERFSX)
void LAPACK_zgerfsx(
    char const* trans, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    double const* R,
    double const* C,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgerq2 LAPACK_GLOBAL(cgerq2,CGERQ2)
void LAPACK_cgerq2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dgerq2 LAPACK_GLOBAL(dgerq2,DGERQ2)
void LAPACK_dgerq2(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work,
    lapack_int* info );

#define LAPACK_sgerq2 LAPACK_GLOBAL(sgerq2,SGERQ2)
void LAPACK_sgerq2(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work,
    lapack_int* info );

#define LAPACK_zgerq2 LAPACK_GLOBAL(zgerq2,ZGERQ2)
void LAPACK_zgerq2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cgerqf LAPACK_GLOBAL(cgerqf,CGERQF)
void LAPACK_cgerqf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgerqf LAPACK_GLOBAL(dgerqf,DGERQF)
void LAPACK_dgerqf(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgerqf LAPACK_GLOBAL(sgerqf,SGERQF)
void LAPACK_sgerqf(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgerqf LAPACK_GLOBAL(zgerqf,ZGERQF)
void LAPACK_zgerqf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgesdd LAPACK_GLOBAL(cgesdd,CGESDD)
void LAPACK_cgesdd(
    char const* jobz,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* S,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* VT, lapack_int const* ldvt,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dgesdd LAPACK_GLOBAL(dgesdd,DGESDD)
void LAPACK_dgesdd(
    char const* jobz,
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* S,
    double* U, lapack_int const* ldu,
    double* VT, lapack_int const* ldvt,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgesdd LAPACK_GLOBAL(sgesdd,SGESDD)
void LAPACK_sgesdd(
    char const* jobz,
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* S,
    float* U, lapack_int const* ldu,
    float* VT, lapack_int const* ldvt,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgesdd LAPACK_GLOBAL(zgesdd,ZGESDD)
void LAPACK_zgesdd(
    char const* jobz,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* S,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* VT, lapack_int const* ldvt,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cgesv LAPACK_GLOBAL(cgesv,CGESV)
void LAPACK_cgesv(
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dgesv LAPACK_GLOBAL(dgesv,DGESV)
void LAPACK_dgesv(
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sgesv LAPACK_GLOBAL(sgesv,SGESV)
void LAPACK_sgesv(
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zgesv LAPACK_GLOBAL(zgesv,ZGESV)
void LAPACK_zgesv(
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsgesv LAPACK_GLOBAL(dsgesv,DSGESV)
void LAPACK_dsgesv(
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* work,
    float* swork, lapack_int* iter,
    lapack_int* info );

#define LAPACK_zcgesv LAPACK_GLOBAL(zcgesv,ZCGESV)
void LAPACK_zcgesv(
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    lapack_complex_double* work,
    lapack_complex_float* swork,
    double* rwork, lapack_int* iter,
    lapack_int* info );

#define LAPACK_cgesvd LAPACK_GLOBAL(cgesvd,CGESVD)
void LAPACK_cgesvd(
    char const* jobu, char const* jobvt,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* S,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* VT, lapack_int const* ldvt,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgesvd LAPACK_GLOBAL(dgesvd,DGESVD)
void LAPACK_dgesvd(
    char const* jobu, char const* jobvt,
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* S,
    double* U, lapack_int const* ldu,
    double* VT, lapack_int const* ldvt,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgesvd LAPACK_GLOBAL(sgesvd,SGESVD)
void LAPACK_sgesvd(
    char const* jobu, char const* jobvt,
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* S,
    float* U, lapack_int const* ldu,
    float* VT, lapack_int const* ldvt,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgesvd LAPACK_GLOBAL(zgesvd,ZGESVD)
void LAPACK_zgesvd(
    char const* jobu, char const* jobvt,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* S,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* VT, lapack_int const* ldvt,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgesvdq LAPACK_GLOBAL(cgesvdq,CGESVDQ)
void LAPACK_cgesvdq(
    char const* joba, char const* jobp, char const* jobr, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* S,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* V, lapack_int const* ldv, lapack_int* numrank,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_complex_float* cwork, lapack_int* lcwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_dgesvdq LAPACK_GLOBAL(dgesvdq,DGESVDQ)
void LAPACK_dgesvdq(
    char const* joba, char const* jobp, char const* jobr, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* S,
    double* U, lapack_int const* ldu,
    double* V, lapack_int const* ldv, lapack_int* numrank,
    lapack_int* iwork, lapack_int const* liwork,
    double* work, lapack_int* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_sgesvdq LAPACK_GLOBAL(sgesvdq,SGESVDQ)
void LAPACK_sgesvdq(
    char const* joba, char const* jobp, char const* jobr, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* S,
    float* U, lapack_int const* ldu,
    float* V, lapack_int const* ldv, lapack_int* numrank,
    lapack_int* iwork, lapack_int const* liwork,
    float* work, lapack_int* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_zgesvdq LAPACK_GLOBAL(zgesvdq,ZGESVDQ)
void LAPACK_zgesvdq(
    char const* joba, char const* jobp, char const* jobr, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* S,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* V, lapack_int const* ldv, lapack_int* numrank,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_complex_double* cwork, lapack_int* lcwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_cgesvdx LAPACK_GLOBAL(cgesvdx,CGESVDX)
void LAPACK_cgesvdx(
    char const* jobu, char const* jobvt, char const* range,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* ns,
    float* S,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* VT, lapack_int const* ldvt,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dgesvdx LAPACK_GLOBAL(dgesvdx,DGESVDX)
void LAPACK_dgesvdx(
    char const* jobu, char const* jobvt, char const* range,
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* ns,
    double* S,
    double* U, lapack_int const* ldu,
    double* VT, lapack_int const* ldvt,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgesvdx LAPACK_GLOBAL(sgesvdx,SGESVDX)
void LAPACK_sgesvdx(
    char const* jobu, char const* jobvt, char const* range,
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* ns,
    float* S,
    float* U, lapack_int const* ldu,
    float* VT, lapack_int const* ldvt,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgesvdx LAPACK_GLOBAL(zgesvdx,ZGESVDX)
void LAPACK_zgesvdx(
    char const* jobu, char const* jobvt, char const* range,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* ns,
    double* S,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* VT, lapack_int const* ldvt,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cgesvj LAPACK_GLOBAL(cgesvj,CGESVJ)
void LAPACK_cgesvj(
    char const* joba, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* SVA, lapack_int const* mv,
    lapack_complex_float* V, lapack_int const* ldv,
    lapack_complex_float* cwork, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_dgesvj LAPACK_GLOBAL(dgesvj,DGESVJ)
void LAPACK_dgesvj(
    char const* joba, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* SVA, lapack_int const* mv,
    double* V, lapack_int const* ldv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgesvj LAPACK_GLOBAL(sgesvj,SGESVJ)
void LAPACK_sgesvj(
    char const* joba, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* SVA, lapack_int const* mv,
    float* V, lapack_int const* ldv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgesvj LAPACK_GLOBAL(zgesvj,ZGESVJ)
void LAPACK_zgesvj(
    char const* joba, char const* jobu, char const* jobv,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* SVA, lapack_int const* mv,
    lapack_complex_double* V, lapack_int const* ldv,
    lapack_complex_double* cwork, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_cgesvx LAPACK_GLOBAL(cgesvx,CGESVX)
void LAPACK_cgesvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgesvx LAPACK_GLOBAL(dgesvx,DGESVX)
void LAPACK_dgesvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgesvx LAPACK_GLOBAL(sgesvx,SGESVX)
void LAPACK_sgesvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgesvx LAPACK_GLOBAL(zgesvx,ZGESVX)
void LAPACK_zgesvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgesvxx LAPACK_GLOBAL(cgesvxx,CGESVXX)
void LAPACK_cgesvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgesvxx LAPACK_GLOBAL(dgesvxx,DGESVXX)
void LAPACK_dgesvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgesvxx LAPACK_GLOBAL(sgesvxx,SGESVXX)
void LAPACK_sgesvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    float* R,
    float* C,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgesvxx LAPACK_GLOBAL(zgesvxx,ZGESVXX)
void LAPACK_zgesvxx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    double* R,
    double* C,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgetf2 LAPACK_GLOBAL(cgetf2,CGETF2)
void LAPACK_cgetf2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_dgetf2 LAPACK_GLOBAL(dgetf2,DGETF2)
void LAPACK_dgetf2(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_sgetf2 LAPACK_GLOBAL(sgetf2,SGETF2)
void LAPACK_sgetf2(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_zgetf2 LAPACK_GLOBAL(zgetf2,ZGETF2)
void LAPACK_zgetf2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_cgetrf LAPACK_GLOBAL(cgetrf,CGETRF)
void LAPACK_cgetrf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_dgetrf LAPACK_GLOBAL(dgetrf,DGETRF)
void LAPACK_dgetrf(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_sgetrf LAPACK_GLOBAL(sgetrf,SGETRF)
void LAPACK_sgetrf(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_zgetrf LAPACK_GLOBAL(zgetrf,ZGETRF)
void LAPACK_zgetrf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_cgetrf2 LAPACK_GLOBAL(cgetrf2,CGETRF2)
void LAPACK_cgetrf2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_dgetrf2 LAPACK_GLOBAL(dgetrf2,DGETRF2)
void LAPACK_dgetrf2(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_sgetrf2 LAPACK_GLOBAL(sgetrf2,SGETRF2)
void LAPACK_sgetrf2(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_zgetrf2 LAPACK_GLOBAL(zgetrf2,ZGETRF2)
void LAPACK_zgetrf2(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_cgetri LAPACK_GLOBAL(cgetri,CGETRI)
void LAPACK_cgetri(
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgetri LAPACK_GLOBAL(dgetri,DGETRI)
void LAPACK_dgetri(
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int const* ipiv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgetri LAPACK_GLOBAL(sgetri,SGETRI)
void LAPACK_sgetri(
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int const* ipiv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgetri LAPACK_GLOBAL(zgetri,ZGETRI)
void LAPACK_zgetri(
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgetrs LAPACK_GLOBAL(cgetrs,CGETRS)
void LAPACK_cgetrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dgetrs LAPACK_GLOBAL(dgetrs,DGETRS)
void LAPACK_dgetrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sgetrs LAPACK_GLOBAL(sgetrs,SGETRS)
void LAPACK_sgetrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zgetrs LAPACK_GLOBAL(zgetrs,ZGETRS)
void LAPACK_zgetrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cgetsls LAPACK_GLOBAL(cgetsls,CGETSLS)
void LAPACK_cgetsls(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgetsls LAPACK_GLOBAL(dgetsls,DGETSLS)
void LAPACK_dgetsls(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgetsls LAPACK_GLOBAL(sgetsls,SGETSLS)
void LAPACK_sgetsls(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgetsls LAPACK_GLOBAL(zgetsls,ZGETSLS)
void LAPACK_zgetsls(
    char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cggbak LAPACK_GLOBAL(cggbak,CGGBAK)
void LAPACK_cggbak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float const* lscale,
    float const* rscale, lapack_int const* m,
    lapack_complex_float* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_dggbak LAPACK_GLOBAL(dggbak,DGGBAK)
void LAPACK_dggbak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double const* lscale,
    double const* rscale, lapack_int const* m,
    double* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_sggbak LAPACK_GLOBAL(sggbak,SGGBAK)
void LAPACK_sggbak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float const* lscale,
    float const* rscale, lapack_int const* m,
    float* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_zggbak LAPACK_GLOBAL(zggbak,ZGGBAK)
void LAPACK_zggbak(
    char const* job, char const* side,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double const* lscale,
    double const* rscale, lapack_int const* m,
    lapack_complex_double* V, lapack_int const* ldv,
    lapack_int* info );

#define LAPACK_cggbal LAPACK_GLOBAL(cggbal,CGGBAL)
void LAPACK_cggbal(
    char const* job,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb, lapack_int* ilo, lapack_int* ihi,
    float* lscale,
    float* rscale,
    float* work,
    lapack_int* info );

#define LAPACK_dggbal LAPACK_GLOBAL(dggbal,DGGBAL)
void LAPACK_dggbal(
    char const* job,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb, lapack_int* ilo, lapack_int* ihi,
    double* lscale,
    double* rscale,
    double* work,
    lapack_int* info );

#define LAPACK_sggbal LAPACK_GLOBAL(sggbal,SGGBAL)
void LAPACK_sggbal(
    char const* job,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb, lapack_int* ilo, lapack_int* ihi,
    float* lscale,
    float* rscale,
    float* work,
    lapack_int* info );

#define LAPACK_zggbal LAPACK_GLOBAL(zggbal,ZGGBAL)
void LAPACK_zggbal(
    char const* job,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb, lapack_int* ilo, lapack_int* ihi,
    double* lscale,
    double* rscale,
    double* work,
    lapack_int* info );

#define LAPACK_cgges LAPACK_GLOBAL(cgges,CGGES)
void LAPACK_cgges(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_C_SELECT2 selctg,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb, lapack_int* sdim,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* VSL, lapack_int const* ldvsl,
    lapack_complex_float* VSR, lapack_int const* ldvsr,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_dgges LAPACK_GLOBAL(dgges,DGGES)
void LAPACK_dgges(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_D_SELECT3 selctg,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb, lapack_int* sdim,
    double* alphar,
    double* alphai,
    double* beta,
    double* VSL, lapack_int const* ldvsl,
    double* VSR, lapack_int const* ldvsr,
    double* work, lapack_int const* lwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_sgges LAPACK_GLOBAL(sgges,SGGES)
void LAPACK_sgges(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_S_SELECT3 selctg,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb, lapack_int* sdim,
    float* alphar,
    float* alphai,
    float* beta,
    float* VSL, lapack_int const* ldvsl,
    float* VSR, lapack_int const* ldvsr,
    float* work, lapack_int const* lwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_zgges LAPACK_GLOBAL(zgges,ZGGES)
void LAPACK_zgges(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_Z_SELECT2 selctg,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb, lapack_int* sdim,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* VSL, lapack_int const* ldvsl,
    lapack_complex_double* VSR, lapack_int const* ldvsr,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_cgges3 LAPACK_GLOBAL(cgges3,CGGES3)
void LAPACK_cgges3(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_C_SELECT2 selctg,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb, lapack_int* sdim,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* VSL, lapack_int const* ldvsl,
    lapack_complex_float* VSR, lapack_int const* ldvsr,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_dgges3 LAPACK_GLOBAL(dgges3,DGGES3)
void LAPACK_dgges3(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_D_SELECT3 selctg,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb, lapack_int* sdim,
    double* alphar,
    double* alphai,
    double* beta,
    double* VSL, lapack_int const* ldvsl,
    double* VSR, lapack_int const* ldvsr,
    double* work, lapack_int const* lwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_sgges3 LAPACK_GLOBAL(sgges3,SGGES3)
void LAPACK_sgges3(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_S_SELECT3 selctg,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb, lapack_int* sdim,
    float* alphar,
    float* alphai,
    float* beta,
    float* VSL, lapack_int const* ldvsl,
    float* VSR, lapack_int const* ldvsr,
    float* work, lapack_int const* lwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_zgges3 LAPACK_GLOBAL(zgges3,ZGGES3)
void LAPACK_zgges3(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_Z_SELECT2 selctg,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb, lapack_int* sdim,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* VSL, lapack_int const* ldvsl,
    lapack_complex_double* VSR, lapack_int const* ldvsr,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_cggesx LAPACK_GLOBAL(cggesx,CGGESX)
void LAPACK_cggesx(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_C_SELECT2 selctg, char const* sense,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb, lapack_int* sdim,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* VSL, lapack_int const* ldvsl,
    lapack_complex_float* VSR, lapack_int const* ldvsr,
    float* rconde,
    float* rcondv,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork, lapack_int const* liwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_dggesx LAPACK_GLOBAL(dggesx,DGGESX)
void LAPACK_dggesx(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_D_SELECT3 selctg, char const* sense,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb, lapack_int* sdim,
    double* alphar,
    double* alphai,
    double* beta,
    double* VSL, lapack_int const* ldvsl,
    double* VSR, lapack_int const* ldvsr,
    double* rconde,
    double* rcondv,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_sggesx LAPACK_GLOBAL(sggesx,SGGESX)
void LAPACK_sggesx(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_S_SELECT3 selctg, char const* sense,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb, lapack_int* sdim,
    float* alphar,
    float* alphai,
    float* beta,
    float* VSL, lapack_int const* ldvsl,
    float* VSR, lapack_int const* ldvsr,
    float* rconde,
    float* rcondv,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_zggesx LAPACK_GLOBAL(zggesx,ZGGESX)
void LAPACK_zggesx(
    char const* jobvsl, char const* jobvsr, char const* sort, LAPACK_Z_SELECT2 selctg, char const* sense,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb, lapack_int* sdim,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* VSL, lapack_int const* ldvsl,
    lapack_complex_double* VSR, lapack_int const* ldvsr,
    double* rconde,
    double* rcondv,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork, lapack_int const* liwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_cggev LAPACK_GLOBAL(cggev,CGGEV)
void LAPACK_cggev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dggev LAPACK_GLOBAL(dggev,DGGEV)
void LAPACK_dggev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* alphar,
    double* alphai,
    double* beta,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sggev LAPACK_GLOBAL(sggev,SGGEV)
void LAPACK_sggev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* alphar,
    float* alphai,
    float* beta,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zggev LAPACK_GLOBAL(zggev,ZGGEV)
void LAPACK_zggev(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cggev3 LAPACK_GLOBAL(cggev3,CGGEV3)
void LAPACK_cggev3(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dggev3 LAPACK_GLOBAL(dggev3,DGGEV3)
void LAPACK_dggev3(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* alphar,
    double* alphai,
    double* beta,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sggev3 LAPACK_GLOBAL(sggev3,SGGEV3)
void LAPACK_sggev3(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* alphar,
    float* alphai,
    float* beta,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zggev3 LAPACK_GLOBAL(zggev3,ZGGEV3)
void LAPACK_zggev3(
    char const* jobvl, char const* jobvr,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cggevx LAPACK_GLOBAL(cggevx,CGGEVX)
void LAPACK_cggevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    float* lscale,
    float* rscale,
    float* abnrm,
    float* bbnrm,
    float* rconde,
    float* rcondv,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_dggevx LAPACK_GLOBAL(dggevx,DGGEVX)
void LAPACK_dggevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* alphar,
    double* alphai,
    double* beta,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    double* lscale,
    double* rscale,
    double* abnrm,
    double* bbnrm,
    double* rconde,
    double* rcondv,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_sggevx LAPACK_GLOBAL(sggevx,SGGEVX)
void LAPACK_sggevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* alphar,
    float* alphai,
    float* beta,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    float* lscale,
    float* rscale,
    float* abnrm,
    float* bbnrm,
    float* rconde,
    float* rcondv,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_zggevx LAPACK_GLOBAL(zggevx,ZGGEVX)
void LAPACK_zggevx(
    char const* balanc, char const* jobvl, char const* jobvr, char const* sense,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr, lapack_int* ilo, lapack_int* ihi,
    double* lscale,
    double* rscale,
    double* abnrm,
    double* bbnrm,
    double* rconde,
    double* rcondv,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork, lapack_logical* BWORK,
    lapack_int* info );

#define LAPACK_cggglm LAPACK_GLOBAL(cggglm,CGGGLM)
void LAPACK_cggglm(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* D,
    lapack_complex_float* X,
    lapack_complex_float* Y,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dggglm LAPACK_GLOBAL(dggglm,DGGGLM)
void LAPACK_dggglm(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* D,
    double* X,
    double* Y,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sggglm LAPACK_GLOBAL(sggglm,SGGGLM)
void LAPACK_sggglm(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* D,
    float* X,
    float* Y,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zggglm LAPACK_GLOBAL(zggglm,ZGGGLM)
void LAPACK_zggglm(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* D,
    lapack_complex_double* X,
    lapack_complex_double* Y,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgghd3 LAPACK_GLOBAL(cgghd3,CGGHD3)
void LAPACK_cgghd3(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgghd3 LAPACK_GLOBAL(dgghd3,DGGHD3)
void LAPACK_dgghd3(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* Q, lapack_int const* ldq,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgghd3 LAPACK_GLOBAL(sgghd3,SGGHD3)
void LAPACK_sgghd3(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* Q, lapack_int const* ldq,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgghd3 LAPACK_GLOBAL(zgghd3,ZGGHD3)
void LAPACK_zgghd3(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgghrd LAPACK_GLOBAL(cgghrd,CGGHRD)
void LAPACK_cgghrd(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_int* info );

#define LAPACK_dgghrd LAPACK_GLOBAL(dgghrd,DGGHRD)
void LAPACK_dgghrd(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* Q, lapack_int const* ldq,
    double* Z, lapack_int const* ldz,
    lapack_int* info );

#define LAPACK_sgghrd LAPACK_GLOBAL(sgghrd,SGGHRD)
void LAPACK_sgghrd(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* Q, lapack_int const* ldq,
    float* Z, lapack_int const* ldz,
    lapack_int* info );

#define LAPACK_zgghrd LAPACK_GLOBAL(zgghrd,ZGGHRD)
void LAPACK_zgghrd(
    char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_int* info );

#define LAPACK_cgglse LAPACK_GLOBAL(cgglse,CGGLSE)
void LAPACK_cgglse(
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* C,
    lapack_complex_float* D,
    lapack_complex_float* X,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dgglse LAPACK_GLOBAL(dgglse,DGGLSE)
void LAPACK_dgglse(
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* C,
    double* D,
    double* X,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sgglse LAPACK_GLOBAL(sgglse,SGGLSE)
void LAPACK_sgglse(
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* C,
    float* D,
    float* X,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zgglse LAPACK_GLOBAL(zgglse,ZGGLSE)
void LAPACK_zgglse(
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* C,
    lapack_complex_double* D,
    lapack_complex_double* X,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cggqrf LAPACK_GLOBAL(cggqrf,CGGQRF)
void LAPACK_cggqrf(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* taua,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* taub,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dggqrf LAPACK_GLOBAL(dggqrf,DGGQRF)
void LAPACK_dggqrf(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    double* A, lapack_int const* lda,
    double* taua,
    double* B, lapack_int const* ldb,
    double* taub,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sggqrf LAPACK_GLOBAL(sggqrf,SGGQRF)
void LAPACK_sggqrf(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    float* A, lapack_int const* lda,
    float* taua,
    float* B, lapack_int const* ldb,
    float* taub,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zggqrf LAPACK_GLOBAL(zggqrf,ZGGQRF)
void LAPACK_zggqrf(
    lapack_int const* n, lapack_int const* m, lapack_int const* p,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* taua,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* taub,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cggrqf LAPACK_GLOBAL(cggrqf,CGGRQF)
void LAPACK_cggrqf(
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* taua,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* taub,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dggrqf LAPACK_GLOBAL(dggrqf,DGGRQF)
void LAPACK_dggrqf(
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* taua,
    double* B, lapack_int const* ldb,
    double* taub,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sggrqf LAPACK_GLOBAL(sggrqf,SGGRQF)
void LAPACK_sggrqf(
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* taua,
    float* B, lapack_int const* ldb,
    float* taub,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zggrqf LAPACK_GLOBAL(zggrqf,ZGGRQF)
void LAPACK_zggrqf(
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* taua,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* taub,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sggsvd LAPACK_GLOBAL(sggsvd,SGGSVD)
lapack_int LAPACK_sggsvd(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    lapack_int* k, lapack_int* l,
    float* a, lapack_int const* lda,
    float* b, lapack_int const* ldb,
    float* alpha, float* beta,
    float* u, lapack_int const* ldu,
    float* v, lapack_int const* ldv,
    float* q, lapack_int const* ldq,
    float* work, lapack_int* iwork, lapack_int* info );

#define LAPACK_dggsvd LAPACK_GLOBAL(dggsvd,DGGSVD)
lapack_int LAPACK_dggsvd(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    lapack_int* k, lapack_int* l,
    double* a, lapack_int const* lda,
    double* b, lapack_int const* ldb,
    double* alpha, double* beta,
    double* u, lapack_int const* ldu,
    double* v, lapack_int const* ldv,
    double* q, lapack_int const* ldq,
    double* work, lapack_int* iwork, lapack_int* info );

#define LAPACK_cggsvd LAPACK_GLOBAL(cggsvd,CGGSVD)
lapack_int LAPACK_cggsvd(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    lapack_int* k, lapack_int* l,
    lapack_complex_float* a, lapack_int const* lda,
    lapack_complex_float* b, lapack_int const* ldb,
    float* alpha, float* beta,
    lapack_complex_float* u, lapack_int const* ldu,
    lapack_complex_float* v, lapack_int const* ldv,
    lapack_complex_float* q, lapack_int const* ldq,
    lapack_complex_float* work, float* rwork,
    lapack_int* iwork, lapack_int* info );

#define LAPACK_zggsvd LAPACK_GLOBAL(zggsvd,ZGGSVD)
lapack_int LAPACK_zggsvd(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p,
    lapack_int* k, lapack_int* l,
    lapack_complex_double* a, lapack_int const* lda,
    lapack_complex_double* b, lapack_int const* ldb,
    double* alpha, double* beta,
    lapack_complex_double* u, lapack_int const* ldu,
    lapack_complex_double* v, lapack_int const* ldv,
    lapack_complex_double* q, lapack_int const* ldq,
    lapack_complex_double* work, double* rwork,
    lapack_int* iwork, lapack_int* info );

#define LAPACK_cggsvd3 LAPACK_GLOBAL(cggsvd3,CGGSVD3)
void LAPACK_cggsvd3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p, lapack_int* k, lapack_int* l,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float* alpha,
    float* beta,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* V, lapack_int const* ldv,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dggsvd3 LAPACK_GLOBAL(dggsvd3,DGGSVD3)
void LAPACK_dggsvd3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p, lapack_int* k, lapack_int* l,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* alpha,
    double* beta,
    double* U, lapack_int const* ldu,
    double* V, lapack_int const* ldv,
    double* Q, lapack_int const* ldq,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sggsvd3 LAPACK_GLOBAL(sggsvd3,SGGSVD3)
void LAPACK_sggsvd3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p, lapack_int* k, lapack_int* l,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* alpha,
    float* beta,
    float* U, lapack_int const* ldu,
    float* V, lapack_int const* ldv,
    float* Q, lapack_int const* ldq,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zggsvd3 LAPACK_GLOBAL(zggsvd3,ZGGSVD3)
void LAPACK_zggsvd3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* n, lapack_int const* p, lapack_int* k, lapack_int* l,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double* alpha,
    double* beta,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* V, lapack_int const* ldv,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sggsvp LAPACK_GLOBAL(sggsvp,SGGSVP)
lapack_int LAPACK_sggsvp(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    float* a, lapack_int const* lda,
    float* b, lapack_int const* ldb,
    float* tola, float* tolb,
    lapack_int* k, lapack_int* l,
    float* u, lapack_int const* ldu,
    float* v, lapack_int const* ldv,
    float* q, lapack_int const* ldq,
    lapack_int* iwork, float* tau,
    float* work, lapack_int* info );

#define LAPACK_dggsvp LAPACK_GLOBAL(dggsvp,DGGSVP)
lapack_int LAPACK_dggsvp(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    double* a, lapack_int const* lda,
    double* b, lapack_int const* ldb,
    double* tola, double* tolb,
    lapack_int* k, lapack_int* l,
    double* u, lapack_int const* ldu,
    double* v, lapack_int const* ldv,
    double* q, lapack_int const* ldq,
    lapack_int* iwork, double* tau,
    double* work, lapack_int* info );

#define LAPACK_cggsvp LAPACK_GLOBAL(cggsvp,CGGSVP)
lapack_int LAPACK_cggsvp(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    lapack_complex_float* a, lapack_int const* lda,
    lapack_complex_float* b, lapack_int const* ldb,
    float* tola, float* tolb, lapack_int* k, lapack_int* l,
    lapack_complex_float* u, lapack_int const* ldu,
    lapack_complex_float* v, lapack_int const* ldv,
    lapack_complex_float* q, lapack_int const* ldq,
    lapack_int* iwork, float* rwork, lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int* info );

#define LAPACK_zggsvp LAPACK_GLOBAL(zggsvp,ZGGSVP)
lapack_int LAPACK_zggsvp(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    lapack_complex_double* a, lapack_int const* lda,
    lapack_complex_double* b, lapack_int const* ldb,
    double* tola, double* tolb, lapack_int* k, lapack_int* l,
    lapack_complex_double* u, lapack_int const* ldu,
    lapack_complex_double* v, lapack_int const* ldv,
    lapack_complex_double* q, lapack_int const* ldq,
    lapack_int* iwork, double* rwork, lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int* info );

#define LAPACK_cggsvp3 LAPACK_GLOBAL(cggsvp3,CGGSVP3)
void LAPACK_cggsvp3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float const* tola,
    float const* tolb, lapack_int* k, lapack_int* l,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* V, lapack_int const* ldv,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_int* iwork,
    float* rwork,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dggsvp3 LAPACK_GLOBAL(dggsvp3,DGGSVP3)
void LAPACK_dggsvp3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double const* tola,
    double const* tolb, lapack_int* k, lapack_int* l,
    double* U, lapack_int const* ldu,
    double* V, lapack_int const* ldv,
    double* Q, lapack_int const* ldq,
    lapack_int* iwork,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sggsvp3 LAPACK_GLOBAL(sggsvp3,SGGSVP3)
void LAPACK_sggsvp3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float const* tola,
    float const* tolb, lapack_int* k, lapack_int* l,
    float* U, lapack_int const* ldu,
    float* V, lapack_int const* ldv,
    float* Q, lapack_int const* ldq,
    lapack_int* iwork,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zggsvp3 LAPACK_GLOBAL(zggsvp3,ZGGSVP3)
void LAPACK_zggsvp3(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double const* tola,
    double const* tolb, lapack_int* k, lapack_int* l,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* V, lapack_int const* ldv,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_int* iwork,
    double* rwork,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cgtcon LAPACK_GLOBAL(cgtcon,CGTCON)
void LAPACK_cgtcon(
    char const* norm,
    lapack_int const* n,
    lapack_complex_float const* DL,
    lapack_complex_float const* D,
    lapack_complex_float const* DU,
    lapack_complex_float const* DU2, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dgtcon LAPACK_GLOBAL(dgtcon,DGTCON)
void LAPACK_dgtcon(
    char const* norm,
    lapack_int const* n,
    double const* DL,
    double const* D,
    double const* DU,
    double const* DU2, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgtcon LAPACK_GLOBAL(sgtcon,SGTCON)
void LAPACK_sgtcon(
    char const* norm,
    lapack_int const* n,
    float const* DL,
    float const* D,
    float const* DU,
    float const* DU2, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgtcon LAPACK_GLOBAL(zgtcon,ZGTCON)
void LAPACK_zgtcon(
    char const* norm,
    lapack_int const* n,
    lapack_complex_double const* DL,
    lapack_complex_double const* D,
    lapack_complex_double const* DU,
    lapack_complex_double const* DU2, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cgtrfs LAPACK_GLOBAL(cgtrfs,CGTRFS)
void LAPACK_cgtrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* DL,
    lapack_complex_float const* D,
    lapack_complex_float const* DU,
    lapack_complex_float const* DLF,
    lapack_complex_float const* DF,
    lapack_complex_float const* DUF,
    lapack_complex_float const* DU2, lapack_int const* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgtrfs LAPACK_GLOBAL(dgtrfs,DGTRFS)
void LAPACK_dgtrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    double const* DL,
    double const* D,
    double const* DU,
    double const* DLF,
    double const* DF,
    double const* DUF,
    double const* DU2, lapack_int const* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgtrfs LAPACK_GLOBAL(sgtrfs,SGTRFS)
void LAPACK_sgtrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    float const* DL,
    float const* D,
    float const* DU,
    float const* DLF,
    float const* DF,
    float const* DUF,
    float const* DU2, lapack_int const* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgtrfs LAPACK_GLOBAL(zgtrfs,ZGTRFS)
void LAPACK_zgtrfs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* DL,
    lapack_complex_double const* D,
    lapack_complex_double const* DU,
    lapack_complex_double const* DLF,
    lapack_complex_double const* DF,
    lapack_complex_double const* DUF,
    lapack_complex_double const* DU2, lapack_int const* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgtsv LAPACK_GLOBAL(cgtsv,CGTSV)
void LAPACK_cgtsv(
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* DL,
    lapack_complex_float* D,
    lapack_complex_float* DU,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dgtsv LAPACK_GLOBAL(dgtsv,DGTSV)
void LAPACK_dgtsv(
    lapack_int const* n, lapack_int const* nrhs,
    double* DL,
    double* D,
    double* DU,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sgtsv LAPACK_GLOBAL(sgtsv,SGTSV)
void LAPACK_sgtsv(
    lapack_int const* n, lapack_int const* nrhs,
    float* DL,
    float* D,
    float* DU,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zgtsv LAPACK_GLOBAL(zgtsv,ZGTSV)
void LAPACK_zgtsv(
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* DL,
    lapack_complex_double* D,
    lapack_complex_double* DU,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cgtsvx LAPACK_GLOBAL(cgtsvx,CGTSVX)
void LAPACK_cgtsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* DL,
    lapack_complex_float const* D,
    lapack_complex_float const* DU,
    lapack_complex_float* DLF,
    lapack_complex_float* DF,
    lapack_complex_float* DUF,
    lapack_complex_float* DU2, lapack_int* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dgtsvx LAPACK_GLOBAL(dgtsvx,DGTSVX)
void LAPACK_dgtsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    double const* DL,
    double const* D,
    double const* DU,
    double* DLF,
    double* DF,
    double* DUF,
    double* DU2, lapack_int* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sgtsvx LAPACK_GLOBAL(sgtsvx,SGTSVX)
void LAPACK_sgtsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    float const* DL,
    float const* D,
    float const* DU,
    float* DLF,
    float* DF,
    float* DUF,
    float* DU2, lapack_int* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zgtsvx LAPACK_GLOBAL(zgtsvx,ZGTSVX)
void LAPACK_zgtsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* DL,
    lapack_complex_double const* D,
    lapack_complex_double const* DU,
    lapack_complex_double* DLF,
    lapack_complex_double* DF,
    lapack_complex_double* DUF,
    lapack_complex_double* DU2, lapack_int* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cgttrf LAPACK_GLOBAL(cgttrf,CGTTRF)
void LAPACK_cgttrf(
    lapack_int const* n,
    lapack_complex_float* DL,
    lapack_complex_float* D,
    lapack_complex_float* DU,
    lapack_complex_float* DU2, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_dgttrf LAPACK_GLOBAL(dgttrf,DGTTRF)
void LAPACK_dgttrf(
    lapack_int const* n,
    double* DL,
    double* D,
    double* DU,
    double* DU2, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_sgttrf LAPACK_GLOBAL(sgttrf,SGTTRF)
void LAPACK_sgttrf(
    lapack_int const* n,
    float* DL,
    float* D,
    float* DU,
    float* DU2, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_zgttrf LAPACK_GLOBAL(zgttrf,ZGTTRF)
void LAPACK_zgttrf(
    lapack_int const* n,
    lapack_complex_double* DL,
    lapack_complex_double* D,
    lapack_complex_double* DU,
    lapack_complex_double* DU2, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_cgttrs LAPACK_GLOBAL(cgttrs,CGTTRS)
void LAPACK_cgttrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* DL,
    lapack_complex_float const* D,
    lapack_complex_float const* DU,
    lapack_complex_float const* DU2, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dgttrs LAPACK_GLOBAL(dgttrs,DGTTRS)
void LAPACK_dgttrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    double const* DL,
    double const* D,
    double const* DU,
    double const* DU2, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sgttrs LAPACK_GLOBAL(sgttrs,SGTTRS)
void LAPACK_sgttrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    float const* DL,
    float const* D,
    float const* DU,
    float const* DU2, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zgttrs LAPACK_GLOBAL(zgttrs,ZGTTRS)
void LAPACK_zgttrs(
    char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* DL,
    lapack_complex_double const* D,
    lapack_complex_double const* DU,
    lapack_complex_double const* DU2, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chbev LAPACK_GLOBAL(chbev,CHBEV)
void LAPACK_chbev(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhbev LAPACK_GLOBAL(zhbev,ZHBEV)
void LAPACK_zhbev(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chbev_2stage LAPACK_GLOBAL(chbev_2stage,CHBEV_2STAGE)
void LAPACK_chbev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhbev_2stage LAPACK_GLOBAL(zhbev_2stage,ZHBEV_2STAGE)
void LAPACK_zhbev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_chbevd LAPACK_GLOBAL(chbevd,CHBEVD)
void LAPACK_chbevd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zhbevd LAPACK_GLOBAL(zhbevd,ZHBEVD)
void LAPACK_zhbevd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_chbevd_2stage LAPACK_GLOBAL(chbevd_2stage,CHBEVD_2STAGE)
void LAPACK_chbevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zhbevd_2stage LAPACK_GLOBAL(zhbevd_2stage,ZHBEVD_2STAGE)
void LAPACK_zhbevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_chbevx LAPACK_GLOBAL(chbevx,CHBEVX)
void LAPACK_chbevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* Q, lapack_int const* ldq,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zhbevx LAPACK_GLOBAL(zhbevx,ZHBEVX)
void LAPACK_zhbevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* Q, lapack_int const* ldq,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_chbevx_2stage LAPACK_GLOBAL(chbevx_2stage,CHBEVX_2STAGE)
void LAPACK_chbevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* Q, lapack_int const* ldq,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zhbevx_2stage LAPACK_GLOBAL(zhbevx_2stage,ZHBEVX_2STAGE)
void LAPACK_zhbevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* Q, lapack_int const* ldq,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_chbgst LAPACK_GLOBAL(chbgst,CHBGST)
void LAPACK_chbgst(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float const* BB, lapack_int const* ldbb,
    lapack_complex_float* X, lapack_int const* ldx,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhbgst LAPACK_GLOBAL(zhbgst,ZHBGST)
void LAPACK_zhbgst(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double const* BB, lapack_int const* ldbb,
    lapack_complex_double* X, lapack_int const* ldx,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chbgv LAPACK_GLOBAL(chbgv,CHBGV)
void LAPACK_chbgv(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* BB, lapack_int const* ldbb,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhbgv LAPACK_GLOBAL(zhbgv,ZHBGV)
void LAPACK_zhbgv(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* BB, lapack_int const* ldbb,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chbgvd LAPACK_GLOBAL(chbgvd,CHBGVD)
void LAPACK_chbgvd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* BB, lapack_int const* ldbb,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zhbgvd LAPACK_GLOBAL(zhbgvd,ZHBGVD)
void LAPACK_zhbgvd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* BB, lapack_int const* ldbb,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_chbgvx LAPACK_GLOBAL(chbgvx,CHBGVX)
void LAPACK_chbgvx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* BB, lapack_int const* ldbb,
    lapack_complex_float* Q, lapack_int const* ldq,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zhbgvx LAPACK_GLOBAL(zhbgvx,ZHBGVX)
void LAPACK_zhbgvx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* BB, lapack_int const* ldbb,
    lapack_complex_double* Q, lapack_int const* ldq,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_chbtrd LAPACK_GLOBAL(chbtrd,CHBTRD)
void LAPACK_chbtrd(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    float* D,
    float* E,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zhbtrd LAPACK_GLOBAL(zhbtrd,ZHBTRD)
void LAPACK_zhbtrd(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    double* D,
    double* E,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_checon LAPACK_GLOBAL(checon,CHECON)
void LAPACK_checon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zhecon LAPACK_GLOBAL(zhecon,ZHECON)
void LAPACK_zhecon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_checon_3 LAPACK_GLOBAL(checon_3,CHECON_3)
void LAPACK_checon_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* E, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zhecon_3 LAPACK_GLOBAL(zhecon_3,ZHECON_3)
void LAPACK_zhecon_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* E, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cheequb LAPACK_GLOBAL(cheequb,CHEEQUB)
void LAPACK_cheequb(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* S,
    float* scond,
    float* amax,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zheequb LAPACK_GLOBAL(zheequb,ZHEEQUB)
void LAPACK_zheequb(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* S,
    double* scond,
    double* amax,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cheev LAPACK_GLOBAL(cheev,CHEEV)
void LAPACK_cheev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* W,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_zheev LAPACK_GLOBAL(zheev,ZHEEV)
void LAPACK_zheev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* W,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cheev_2stage LAPACK_GLOBAL(cheev_2stage,CHEEV_2STAGE)
void LAPACK_cheev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* W,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_zheev_2stage LAPACK_GLOBAL(zheev_2stage,ZHEEV_2STAGE)
void LAPACK_zheev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* W,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_cheevd LAPACK_GLOBAL(cheevd,CHEEVD)
void LAPACK_cheevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* W,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zheevd LAPACK_GLOBAL(zheevd,ZHEEVD)
void LAPACK_zheevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* W,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_cheevd_2stage LAPACK_GLOBAL(cheevd_2stage,CHEEVD_2STAGE)
void LAPACK_cheevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* W,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zheevd_2stage LAPACK_GLOBAL(zheevd_2stage,ZHEEVD_2STAGE)
void LAPACK_zheevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* W,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_cheevr LAPACK_GLOBAL(cheevr,CHEEVR)
void LAPACK_cheevr(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zheevr LAPACK_GLOBAL(zheevr,ZHEEVR)
void LAPACK_zheevr(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_cheevr_2stage LAPACK_GLOBAL(cheevr_2stage,CHEEVR_2STAGE)
void LAPACK_cheevr_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zheevr_2stage LAPACK_GLOBAL(zheevr_2stage,ZHEEVR_2STAGE)
void LAPACK_zheevr_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_cheevx LAPACK_GLOBAL(cheevx,CHEEVX)
void LAPACK_cheevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zheevx LAPACK_GLOBAL(zheevx,ZHEEVX)
void LAPACK_zheevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_cheevx_2stage LAPACK_GLOBAL(cheevx_2stage,CHEEVX_2STAGE)
void LAPACK_cheevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zheevx_2stage LAPACK_GLOBAL(zheevx_2stage,ZHEEVX_2STAGE)
void LAPACK_zheevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_chegst LAPACK_GLOBAL(chegst,CHEGST)
void LAPACK_chegst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zhegst LAPACK_GLOBAL(zhegst,ZHEGST)
void LAPACK_zhegst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chegv LAPACK_GLOBAL(chegv,CHEGV)
void LAPACK_chegv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float* W,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhegv LAPACK_GLOBAL(zhegv,ZHEGV)
void LAPACK_zhegv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double* W,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_chegv_2stage LAPACK_GLOBAL(chegv_2stage,CHEGV_2STAGE)
void LAPACK_chegv_2stage(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float* W,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhegv_2stage LAPACK_GLOBAL(zhegv_2stage,ZHEGV_2STAGE)
void LAPACK_zhegv_2stage(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double* W,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_chegvd LAPACK_GLOBAL(chegvd,CHEGVD)
void LAPACK_chegvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float* W,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zhegvd LAPACK_GLOBAL(zhegvd,ZHEGVD)
void LAPACK_zhegvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double* W,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_chegvx LAPACK_GLOBAL(chegvx,CHEGVX)
void LAPACK_chegvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zhegvx LAPACK_GLOBAL(zhegvx,ZHEGVX)
void LAPACK_zhegvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_cherfs LAPACK_GLOBAL(cherfs,CHERFS)
void LAPACK_cherfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zherfs LAPACK_GLOBAL(zherfs,ZHERFS)
void LAPACK_zherfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cherfsx LAPACK_GLOBAL(cherfsx,CHERFSX)
void LAPACK_cherfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    float* S,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zherfsx LAPACK_GLOBAL(zherfsx,ZHERFSX)
void LAPACK_zherfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    double* S,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chesv LAPACK_GLOBAL(chesv,CHESV)
void LAPACK_chesv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhesv LAPACK_GLOBAL(zhesv,ZHESV)
void LAPACK_zhesv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chesv_aa LAPACK_GLOBAL(chesv_aa,CHESV_AA)
void LAPACK_chesv_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhesv_aa LAPACK_GLOBAL(zhesv_aa,ZHESV_AA)
void LAPACK_zhesv_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chesv_aa_2stage LAPACK_GLOBAL(chesv_aa_2stage,CHESV_AA_2STAGE)
void LAPACK_chesv_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhesv_aa_2stage LAPACK_GLOBAL(zhesv_aa_2stage,ZHESV_AA_2STAGE)
void LAPACK_zhesv_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chesv_rk LAPACK_GLOBAL(chesv_rk,CHESV_RK)
void LAPACK_chesv_rk(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* E, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhesv_rk LAPACK_GLOBAL(zhesv_rk,ZHESV_RK)
void LAPACK_zhesv_rk(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* E, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chesv_rook LAPACK_GLOBAL(chesv_rook,CHESV_ROOK)
void LAPACK_chesv_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhesv_rook LAPACK_GLOBAL(zhesv_rook,ZHESV_ROOK)
void LAPACK_zhesv_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chesvx LAPACK_GLOBAL(chesvx,CHESVX)
void LAPACK_chesvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, lapack_int* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhesvx LAPACK_GLOBAL(zhesvx,ZHESVX)
void LAPACK_zhesvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, lapack_int* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_chesvxx LAPACK_GLOBAL(chesvxx,CHESVXX)
void LAPACK_chesvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    float* S,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhesvxx LAPACK_GLOBAL(zhesvxx,ZHESVXX)
void LAPACK_zhesvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    double* S,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cheswapr LAPACK_GLOBAL(cheswapr,CHESWAPR)
void LAPACK_cheswapr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* i1, lapack_int const* i2 );

#define LAPACK_zheswapr LAPACK_GLOBAL(zheswapr,ZHESWAPR)
void LAPACK_zheswapr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* i1, lapack_int const* i2 );

#define LAPACK_chetrd LAPACK_GLOBAL(chetrd,CHETRD)
void LAPACK_chetrd(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* D,
    float* E,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrd LAPACK_GLOBAL(zhetrd,ZHETRD)
void LAPACK_zhetrd(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* D,
    double* E,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrd_2stage LAPACK_GLOBAL(chetrd_2stage,CHETRD_2STAGE)
void LAPACK_chetrd_2stage(
    char const* vect, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    float* D,
    float* E,
    lapack_complex_float* tau,
    lapack_complex_float* HOUS2, lapack_int const* lhous2,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrd_2stage LAPACK_GLOBAL(zhetrd_2stage,ZHETRD_2STAGE)
void LAPACK_zhetrd_2stage(
    char const* vect, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    double* D,
    double* E,
    lapack_complex_double* tau,
    lapack_complex_double* HOUS2, lapack_int const* lhous2,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrf LAPACK_GLOBAL(chetrf,CHETRF)
void LAPACK_chetrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrf LAPACK_GLOBAL(zhetrf,ZHETRF)
void LAPACK_zhetrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrf_aa LAPACK_GLOBAL(chetrf_aa,CHETRF_AA)
void LAPACK_chetrf_aa(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrf_aa LAPACK_GLOBAL(zhetrf_aa,ZHETRF_AA)
void LAPACK_zhetrf_aa(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrf_aa_2stage LAPACK_GLOBAL(chetrf_aa_2stage,CHETRF_AA_2STAGE)
void LAPACK_chetrf_aa_2stage(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrf_aa_2stage LAPACK_GLOBAL(zhetrf_aa_2stage,ZHETRF_AA_2STAGE)
void LAPACK_zhetrf_aa_2stage(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrf_rk LAPACK_GLOBAL(chetrf_rk,CHETRF_RK)
void LAPACK_chetrf_rk(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* E, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrf_rk LAPACK_GLOBAL(zhetrf_rk,ZHETRF_RK)
void LAPACK_zhetrf_rk(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* E, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrf_rook LAPACK_GLOBAL(chetrf_rook,CHETRF_ROOK)
void LAPACK_chetrf_rook(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrf_rook LAPACK_GLOBAL(zhetrf_rook,ZHETRF_ROOK)
void LAPACK_zhetrf_rook(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetri LAPACK_GLOBAL(chetri,CHETRI)
void LAPACK_chetri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zhetri LAPACK_GLOBAL(zhetri,ZHETRI)
void LAPACK_zhetri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_chetri2 LAPACK_GLOBAL(chetri2,CHETRI2)
void LAPACK_chetri2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetri2 LAPACK_GLOBAL(zhetri2,ZHETRI2)
void LAPACK_zhetri2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetri2x LAPACK_GLOBAL(chetri2x,CHETRI2X)
void LAPACK_chetri2x(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* work, lapack_int const* nb,
    lapack_int* info );

#define LAPACK_zhetri2x LAPACK_GLOBAL(zhetri2x,ZHETRI2X)
void LAPACK_zhetri2x(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* work, lapack_int const* nb,
    lapack_int* info );

#define LAPACK_chetri_3 LAPACK_GLOBAL(chetri_3,CHETRI_3)
void LAPACK_chetri_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* E, lapack_int const* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetri_3 LAPACK_GLOBAL(zhetri_3,ZHETRI_3)
void LAPACK_zhetri_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* E, lapack_int const* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrs LAPACK_GLOBAL(chetrs,CHETRS)
void LAPACK_chetrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zhetrs LAPACK_GLOBAL(zhetrs,ZHETRS)
void LAPACK_zhetrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chetrs2 LAPACK_GLOBAL(chetrs2,CHETRS2)
void LAPACK_chetrs2(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zhetrs2 LAPACK_GLOBAL(zhetrs2,ZHETRS2)
void LAPACK_zhetrs2(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_chetrs_3 LAPACK_GLOBAL(chetrs_3,CHETRS_3)
void LAPACK_chetrs_3(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* E, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zhetrs_3 LAPACK_GLOBAL(zhetrs_3,ZHETRS_3)
void LAPACK_zhetrs_3(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* E, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chetrs_aa LAPACK_GLOBAL(chetrs_aa,CHETRS_AA)
void LAPACK_chetrs_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhetrs_aa LAPACK_GLOBAL(zhetrs_aa,ZHETRS_AA)
void LAPACK_zhetrs_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_chetrs_aa_2stage LAPACK_GLOBAL(chetrs_aa_2stage,CHETRS_AA_2STAGE)
void LAPACK_chetrs_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* TB, lapack_int const* ltb, lapack_int const* ipiv, lapack_int const* ipiv2,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zhetrs_aa_2stage LAPACK_GLOBAL(zhetrs_aa_2stage,ZHETRS_AA_2STAGE)
void LAPACK_zhetrs_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* TB, lapack_int const* ltb, lapack_int const* ipiv, lapack_int const* ipiv2,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chetrs_rook LAPACK_GLOBAL(chetrs_rook,CHETRS_ROOK)
void LAPACK_chetrs_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zhetrs_rook LAPACK_GLOBAL(zhetrs_rook,ZHETRS_ROOK)
void LAPACK_zhetrs_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chfrk LAPACK_GLOBAL(chfrk,CHFRK)
void LAPACK_chfrk(
    char const* transr, char const* uplo, char const* trans,
    lapack_int const* n, lapack_int const* k,
    float const* alpha,
    lapack_complex_float const* A, lapack_int const* lda,
    float const* beta,
    lapack_complex_float* C );

#define LAPACK_zhfrk LAPACK_GLOBAL(zhfrk,ZHFRK)
void LAPACK_zhfrk(
    char const* transr, char const* uplo, char const* trans,
    lapack_int const* n, lapack_int const* k,
    double const* alpha,
    lapack_complex_double const* A, lapack_int const* lda,
    double const* beta,
    lapack_complex_double* C );

#define LAPACK_chgeqz LAPACK_GLOBAL(chgeqz,CHGEQZ)
void LAPACK_chgeqz(
    char const* job, char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_float* H, lapack_int const* ldh,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dhgeqz LAPACK_GLOBAL(dhgeqz,DHGEQZ)
void LAPACK_dhgeqz(
    char const* job, char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double* H, lapack_int const* ldh,
    double* T, lapack_int const* ldt,
    double* alphar,
    double* alphai,
    double* beta,
    double* Q, lapack_int const* ldq,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_shgeqz LAPACK_GLOBAL(shgeqz,SHGEQZ)
void LAPACK_shgeqz(
    char const* job, char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float* H, lapack_int const* ldh,
    float* T, lapack_int const* ldt,
    float* alphar,
    float* alphai,
    float* beta,
    float* Q, lapack_int const* ldq,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhgeqz LAPACK_GLOBAL(zhgeqz,ZHGEQZ)
void LAPACK_zhgeqz(
    char const* job, char const* compq, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_double* H, lapack_int const* ldh,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_chpcon LAPACK_GLOBAL(chpcon,CHPCON)
void LAPACK_chpcon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zhpcon LAPACK_GLOBAL(zhpcon,ZHPCON)
void LAPACK_zhpcon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_chpev LAPACK_GLOBAL(chpev,CHPEV)
void LAPACK_chpev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhpev LAPACK_GLOBAL(zhpev,ZHPEV)
void LAPACK_zhpev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chpevd LAPACK_GLOBAL(chpevd,CHPEVD)
void LAPACK_chpevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zhpevd LAPACK_GLOBAL(zhpevd,ZHPEVD)
void LAPACK_zhpevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_chpevx LAPACK_GLOBAL(chpevx,CHPEVX)
void LAPACK_chpevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zhpevx LAPACK_GLOBAL(zhpevx,ZHPEVX)
void LAPACK_zhpevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_chpgst LAPACK_GLOBAL(chpgst,CHPGST)
void LAPACK_chpgst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    lapack_complex_float const* BP,
    lapack_int* info );

#define LAPACK_zhpgst LAPACK_GLOBAL(zhpgst,ZHPGST)
void LAPACK_zhpgst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    lapack_complex_double const* BP,
    lapack_int* info );

#define LAPACK_chpgv LAPACK_GLOBAL(chpgv,CHPGV)
void LAPACK_chpgv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    lapack_complex_float* BP,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhpgv LAPACK_GLOBAL(zhpgv,ZHPGV)
void LAPACK_zhpgv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    lapack_complex_double* BP,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chpgvd LAPACK_GLOBAL(chpgvd,CHPGVD)
void LAPACK_chpgvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    lapack_complex_float* BP,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zhpgvd LAPACK_GLOBAL(zhpgvd,ZHPGVD)
void LAPACK_zhpgvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    lapack_complex_double* BP,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_chpgvx LAPACK_GLOBAL(chpgvx,CHPGVX)
void LAPACK_chpgvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    lapack_complex_float* BP,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zhpgvx LAPACK_GLOBAL(zhpgvx,ZHPGVX)
void LAPACK_zhpgvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    lapack_complex_double* BP,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_chprfs LAPACK_GLOBAL(chprfs,CHPRFS)
void LAPACK_chprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float const* AFP, lapack_int const* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhprfs LAPACK_GLOBAL(zhprfs,ZHPRFS)
void LAPACK_zhprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double const* AFP, lapack_int const* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chpsv LAPACK_GLOBAL(chpsv,CHPSV)
void LAPACK_chpsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* AP, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zhpsv LAPACK_GLOBAL(zhpsv,ZHPSV)
void LAPACK_zhpsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* AP, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chpsvx LAPACK_GLOBAL(chpsvx,CHPSVX)
void LAPACK_chpsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float* AFP, lapack_int* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_zhpsvx LAPACK_GLOBAL(zhpsvx,ZHPSVX)
void LAPACK_zhpsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double* AFP, lapack_int* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_chptrd LAPACK_GLOBAL(chptrd,CHPTRD)
void LAPACK_chptrd(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    float* D,
    float* E,
    lapack_complex_float* tau,
    lapack_int* info );

#define LAPACK_zhptrd LAPACK_GLOBAL(zhptrd,ZHPTRD)
void LAPACK_zhptrd(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    double* D,
    double* E,
    lapack_complex_double* tau,
    lapack_int* info );

#define LAPACK_chptrf LAPACK_GLOBAL(chptrf,CHPTRF)
void LAPACK_chptrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_zhptrf LAPACK_GLOBAL(zhptrf,ZHPTRF)
void LAPACK_zhptrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_chptri LAPACK_GLOBAL(chptri,CHPTRI)
void LAPACK_chptri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP, lapack_int const* ipiv,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zhptri LAPACK_GLOBAL(zhptri,ZHPTRI)
void LAPACK_zhptri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP, lapack_int const* ipiv,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_chptrs LAPACK_GLOBAL(chptrs,CHPTRS)
void LAPACK_chptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zhptrs LAPACK_GLOBAL(zhptrs,ZHPTRS)
void LAPACK_zhptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_chsein LAPACK_GLOBAL(chsein,CHSEIN)
void LAPACK_chsein(
    char const* side, char const* eigsrc, char const* initv,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_float const* H, lapack_int const* ldh,
    lapack_complex_float* W,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_float* work,
    float* rwork, lapack_int* IFAILL, lapack_int* IFAILR,
    lapack_int* info );

#define LAPACK_dhsein LAPACK_GLOBAL(dhsein,DHSEIN)
void LAPACK_dhsein(
    char const* side, char const* eigsrc, char const* initv,
    lapack_logical* select,
    lapack_int const* n,
    double const* H, lapack_int const* ldh,
    double* WR,
    double const* WI,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    double* work, lapack_int* IFAILL, lapack_int* IFAILR,
    lapack_int* info );

#define LAPACK_shsein LAPACK_GLOBAL(shsein,SHSEIN)
void LAPACK_shsein(
    char const* side, char const* eigsrc, char const* initv,
    lapack_logical* select,
    lapack_int const* n,
    float const* H, lapack_int const* ldh,
    float* WR,
    float const* WI,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    float* work, lapack_int* IFAILL, lapack_int* IFAILR,
    lapack_int* info );

#define LAPACK_zhsein LAPACK_GLOBAL(zhsein,ZHSEIN)
void LAPACK_zhsein(
    char const* side, char const* eigsrc, char const* initv,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_double const* H, lapack_int const* ldh,
    lapack_complex_double* W,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_double* work,
    double* rwork, lapack_int* IFAILL, lapack_int* IFAILR,
    lapack_int* info );

#define LAPACK_chseqr LAPACK_GLOBAL(chseqr,CHSEQR)
void LAPACK_chseqr(
    char const* job, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_float* H, lapack_int const* ldh,
    lapack_complex_float* W,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dhseqr LAPACK_GLOBAL(dhseqr,DHSEQR)
void LAPACK_dhseqr(
    char const* job, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double* H, lapack_int const* ldh,
    double* WR,
    double* WI,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_shseqr LAPACK_GLOBAL(shseqr,SHSEQR)
void LAPACK_shseqr(
    char const* job, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float* H, lapack_int const* ldh,
    float* WR,
    float* WI,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zhseqr LAPACK_GLOBAL(zhseqr,ZHSEQR)
void LAPACK_zhseqr(
    char const* job, char const* compz,
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_double* H, lapack_int const* ldh,
    lapack_complex_double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_clacgv LAPACK_GLOBAL(clacgv,CLACGV)
void LAPACK_clacgv(
    lapack_int const* n,
    lapack_complex_float* X, lapack_int const* incx );

#define LAPACK_zlacgv LAPACK_GLOBAL(zlacgv,ZLACGV)
void LAPACK_zlacgv(
    lapack_int const* n,
    lapack_complex_double* X, lapack_int const* incx );

#define LAPACK_clacn2 LAPACK_GLOBAL(clacn2,CLACN2)
void LAPACK_clacn2(
    lapack_int const* n,
    lapack_complex_float* V,
    lapack_complex_float* X,
    float* est, lapack_int* kase, lapack_int* ISAVE );

#define LAPACK_dlacn2 LAPACK_GLOBAL(dlacn2,DLACN2)
void LAPACK_dlacn2(
    lapack_int const* n,
    double* V,
    double* X, lapack_int* ISGN,
    double* est, lapack_int* kase, lapack_int* ISAVE );

#define LAPACK_slacn2 LAPACK_GLOBAL(slacn2,SLACN2)
void LAPACK_slacn2(
    lapack_int const* n,
    float* V,
    float* X, lapack_int* ISGN,
    float* est, lapack_int* kase, lapack_int* ISAVE );

#define LAPACK_zlacn2 LAPACK_GLOBAL(zlacn2,ZLACN2)
void LAPACK_zlacn2(
    lapack_int const* n,
    lapack_complex_double* V,
    lapack_complex_double* X,
    double* est, lapack_int* kase, lapack_int* ISAVE );

#define LAPACK_clacp2 LAPACK_GLOBAL(clacp2,CLACP2)
void LAPACK_clacp2(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb );

#define LAPACK_zlacp2 LAPACK_GLOBAL(zlacp2,ZLACP2)
void LAPACK_zlacp2(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb );

#define LAPACK_clacpy LAPACK_GLOBAL(clacpy,CLACPY)
void LAPACK_clacpy(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb );

#define LAPACK_dlacpy LAPACK_GLOBAL(dlacpy,DLACPY)
void LAPACK_dlacpy(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* B, lapack_int const* ldb );

#define LAPACK_slacpy LAPACK_GLOBAL(slacpy,SLACPY)
void LAPACK_slacpy(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* B, lapack_int const* ldb );

#define LAPACK_zlacpy LAPACK_GLOBAL(zlacpy,ZLACPY)
void LAPACK_zlacpy(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb );

#define LAPACK_clacrm LAPACK_GLOBAL(clacrm,CLACRM)
void LAPACK_clacrm(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float const* B, lapack_int const* ldb,
    lapack_complex_float* C, lapack_int const* ldc,
    float* rwork );

#define LAPACK_zlacrm LAPACK_GLOBAL(zlacrm,ZLACRM)
void LAPACK_zlacrm(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double const* B, lapack_int const* ldb,
    lapack_complex_double* C, lapack_int const* ldc,
    double* rwork );

#define LAPACK_zlag2c LAPACK_GLOBAL(zlag2c,ZLAG2C)
void LAPACK_zlag2c(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_float* SA, lapack_int const* ldsa,
    lapack_int* info );

#define LAPACK_slag2d LAPACK_GLOBAL(slag2d,SLAG2D)
void LAPACK_slag2d(
    lapack_int const* m, lapack_int const* n,
    float const* SA, lapack_int const* ldsa,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dlag2s LAPACK_GLOBAL(dlag2s,DLAG2S)
void LAPACK_dlag2s(
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    float* SA, lapack_int const* ldsa,
    lapack_int* info );

#define LAPACK_clag2z LAPACK_GLOBAL(clag2z,CLAG2Z)
void LAPACK_clag2z(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* SA, lapack_int const* ldsa,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_clagge LAPACK_GLOBAL(clagge,CLAGGE)
void LAPACK_clagge(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    float const* D,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* iseed,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dlagge LAPACK_GLOBAL(dlagge,DLAGGE)
void LAPACK_dlagge(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    double const* D,
    double* A, lapack_int const* lda, lapack_int* iseed,
    double* work,
    lapack_int* info );

#define LAPACK_slagge LAPACK_GLOBAL(slagge,SLAGGE)
void LAPACK_slagge(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    float const* D,
    float* A, lapack_int const* lda, lapack_int* iseed,
    float* work,
    lapack_int* info );

#define LAPACK_zlagge LAPACK_GLOBAL(zlagge,ZLAGGE)
void LAPACK_zlagge(
    lapack_int const* m, lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    double const* D,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* iseed,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_claghe LAPACK_GLOBAL(claghe,CLAGHE)
void LAPACK_claghe(
    lapack_int const* n, lapack_int const* k,
    float const* D,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* iseed,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zlaghe LAPACK_GLOBAL(zlaghe,ZLAGHE)
void LAPACK_zlaghe(
    lapack_int const* n, lapack_int const* k,
    double const* D,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* iseed,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_clagsy LAPACK_GLOBAL(clagsy,CLAGSY)
void LAPACK_clagsy(
    lapack_int const* n, lapack_int const* k,
    float const* D,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* iseed,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dlagsy LAPACK_GLOBAL(dlagsy,DLAGSY)
void LAPACK_dlagsy(
    lapack_int const* n, lapack_int const* k,
    double const* D,
    double* A, lapack_int const* lda, lapack_int* iseed,
    double* work,
    lapack_int* info );

#define LAPACK_slagsy LAPACK_GLOBAL(slagsy,SLAGSY)
void LAPACK_slagsy(
    lapack_int const* n, lapack_int const* k,
    float const* D,
    float* A, lapack_int const* lda, lapack_int* iseed,
    float* work,
    lapack_int* info );

#define LAPACK_zlagsy LAPACK_GLOBAL(zlagsy,ZLAGSY)
void LAPACK_zlagsy(
    lapack_int const* n, lapack_int const* k,
    double const* D,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* iseed,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_dlamch LAPACK_GLOBAL(dlamch,DLAMCH)
double LAPACK_dlamch(
    char const* cmach );

#define LAPACK_slamch LAPACK_GLOBAL(slamch,SLAMCH)
lapack_float_return LAPACK_slamch(
    char const* cmach );

#define LAPACK_clangb LAPACK_GLOBAL(clangb,CLANGB)
lapack_float_return LAPACK_clangb(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* work );

#define LAPACK_dlangb LAPACK_GLOBAL(dlangb,DLANGB)
double LAPACK_dlangb(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    double const* AB, lapack_int const* ldab,
    double* work );

#define LAPACK_slangb LAPACK_GLOBAL(slangb,SLANGB)
lapack_float_return LAPACK_slangb(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    float const* AB, lapack_int const* ldab,
    float* work );

#define LAPACK_zlangb LAPACK_GLOBAL(zlangb,ZLANGB)
double LAPACK_zlangb(
    char const* norm,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* work );

#define LAPACK_clange LAPACK_GLOBAL(clange,CLANGE)
lapack_float_return LAPACK_clange(
    char const* norm,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_dlange LAPACK_GLOBAL(dlange,DLANGE)
double LAPACK_dlange(
    char const* norm,
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_slange LAPACK_GLOBAL(slange,SLANGE)
lapack_float_return LAPACK_slange(
    char const* norm,
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_zlange LAPACK_GLOBAL(zlange,ZLANGE)
double LAPACK_zlange(
    char const* norm,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_clangt LAPACK_GLOBAL(clangt,CLANGT)
lapack_float_return LAPACK_clangt(
    char const* norm,
    lapack_int const* n,
    lapack_complex_float const* DL,
    lapack_complex_float const* D,
    lapack_complex_float const* DU );

#define LAPACK_dlangt LAPACK_GLOBAL(dlangt,DLANGT)
double LAPACK_dlangt(
    char const* norm,
    lapack_int const* n,
    double const* DL,
    double const* D,
    double const* DU );

#define LAPACK_slangt LAPACK_GLOBAL(slangt,SLANGT)
lapack_float_return LAPACK_slangt(
    char const* norm,
    lapack_int const* n,
    float const* DL,
    float const* D,
    float const* DU );

#define LAPACK_zlangt LAPACK_GLOBAL(zlangt,ZLANGT)
double LAPACK_zlangt(
    char const* norm,
    lapack_int const* n,
    lapack_complex_double const* DL,
    lapack_complex_double const* D,
    lapack_complex_double const* DU );

#define LAPACK_clanhb LAPACK_GLOBAL(clanhb,CLANHB)
lapack_float_return LAPACK_clanhb(
    char const* norm, char const* uplo,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* work );

#define LAPACK_zlanhb LAPACK_GLOBAL(zlanhb,ZLANHB)
double LAPACK_zlanhb(
    char const* norm, char const* uplo,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* work );

#define LAPACK_clanhe LAPACK_GLOBAL(clanhe,CLANHE)
lapack_float_return LAPACK_clanhe(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_zlanhe LAPACK_GLOBAL(zlanhe,ZLANHE)
double LAPACK_zlanhe(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_clanhp LAPACK_GLOBAL(clanhp,CLANHP)
lapack_float_return LAPACK_clanhp(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP,
    float* work );

#define LAPACK_zlanhp LAPACK_GLOBAL(zlanhp,ZLANHP)
double LAPACK_zlanhp(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP,
    double* work );

#define LAPACK_clanhs LAPACK_GLOBAL(clanhs,CLANHS)
lapack_float_return LAPACK_clanhs(
    char const* norm,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_dlanhs LAPACK_GLOBAL(dlanhs,DLANHS)
double LAPACK_dlanhs(
    char const* norm,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_slanhs LAPACK_GLOBAL(slanhs,SLANHS)
lapack_float_return LAPACK_slanhs(
    char const* norm,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_zlanhs LAPACK_GLOBAL(zlanhs,ZLANHS)
double LAPACK_zlanhs(
    char const* norm,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_clanht LAPACK_GLOBAL(clanht,CLANHT)
lapack_float_return LAPACK_clanht(
    char const* norm,
    lapack_int const* n,
    float const* D,
    lapack_complex_float const* E );

#define LAPACK_zlanht LAPACK_GLOBAL(zlanht,ZLANHT)
double LAPACK_zlanht(
    char const* norm,
    lapack_int const* n,
    double const* D,
    lapack_complex_double const* E );

#define LAPACK_clansb LAPACK_GLOBAL(clansb,CLANSB)
lapack_float_return LAPACK_clansb(
    char const* norm, char const* uplo,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* work );

#define LAPACK_dlansb LAPACK_GLOBAL(dlansb,DLANSB)
double LAPACK_dlansb(
    char const* norm, char const* uplo,
    lapack_int const* n, lapack_int const* k,
    double const* AB, lapack_int const* ldab,
    double* work );

#define LAPACK_slansb LAPACK_GLOBAL(slansb,SLANSB)
lapack_float_return LAPACK_slansb(
    char const* norm, char const* uplo,
    lapack_int const* n, lapack_int const* k,
    float const* AB, lapack_int const* ldab,
    float* work );

#define LAPACK_zlansb LAPACK_GLOBAL(zlansb,ZLANSB)
double LAPACK_zlansb(
    char const* norm, char const* uplo,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* work );

#define LAPACK_clansp LAPACK_GLOBAL(clansp,CLANSP)
lapack_float_return LAPACK_clansp(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP,
    float* work );

#define LAPACK_dlansp LAPACK_GLOBAL(dlansp,DLANSP)
double LAPACK_dlansp(
    char const* norm, char const* uplo,
    lapack_int const* n,
    double const* AP,
    double* work );

#define LAPACK_slansp LAPACK_GLOBAL(slansp,SLANSP)
lapack_float_return LAPACK_slansp(
    char const* norm, char const* uplo,
    lapack_int const* n,
    float const* AP,
    float* work );

#define LAPACK_zlansp LAPACK_GLOBAL(zlansp,ZLANSP)
double LAPACK_zlansp(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP,
    double* work );

#define LAPACK_dlanst LAPACK_GLOBAL(dlanst,DLANST)
double LAPACK_dlanst(
    char const* norm,
    lapack_int const* n,
    double const* D,
    double const* E );

#define LAPACK_slanst LAPACK_GLOBAL(slanst,SLANST)
lapack_float_return LAPACK_slanst(
    char const* norm,
    lapack_int const* n,
    float const* D,
    float const* E );

#define LAPACK_clansy LAPACK_GLOBAL(clansy,CLANSY)
lapack_float_return LAPACK_clansy(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_dlansy LAPACK_GLOBAL(dlansy,DLANSY)
double LAPACK_dlansy(
    char const* norm, char const* uplo,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_slansy LAPACK_GLOBAL(slansy,SLANSY)
lapack_float_return LAPACK_slansy(
    char const* norm, char const* uplo,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_zlansy LAPACK_GLOBAL(zlansy,ZLANSY)
double LAPACK_zlansy(
    char const* norm, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_clantb LAPACK_GLOBAL(clantb,CLANTB)
lapack_float_return LAPACK_clantb(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* work );

#define LAPACK_dlantb LAPACK_GLOBAL(dlantb,DLANTB)
double LAPACK_dlantb(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* k,
    double const* AB, lapack_int const* ldab,
    double* work );

#define LAPACK_slantb LAPACK_GLOBAL(slantb,SLANTB)
lapack_float_return LAPACK_slantb(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* k,
    float const* AB, lapack_int const* ldab,
    float* work );

#define LAPACK_zlantb LAPACK_GLOBAL(zlantb,ZLANTB)
double LAPACK_zlantb(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* work );

#define LAPACK_clantp LAPACK_GLOBAL(clantp,CLANTP)
lapack_float_return LAPACK_clantp(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_float const* AP,
    float* work );

#define LAPACK_dlantp LAPACK_GLOBAL(dlantp,DLANTP)
double LAPACK_dlantp(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    double const* AP,
    double* work );

#define LAPACK_slantp LAPACK_GLOBAL(slantp,SLANTP)
lapack_float_return LAPACK_slantp(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    float const* AP,
    float* work );

#define LAPACK_zlantp LAPACK_GLOBAL(zlantp,ZLANTP)
double LAPACK_zlantp(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_double const* AP,
    double* work );

#define LAPACK_clantr LAPACK_GLOBAL(clantr,CLANTR)
lapack_float_return LAPACK_clantr(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_dlantr LAPACK_GLOBAL(dlantr,DLANTR)
double LAPACK_dlantr(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_slantr LAPACK_GLOBAL(slantr,SLANTR)
lapack_float_return LAPACK_slantr(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* work );

#define LAPACK_zlantr LAPACK_GLOBAL(zlantr,ZLANTR)
double LAPACK_zlantr(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* work );

#define LAPACK_clapmr LAPACK_GLOBAL(clapmr,CLAPMR)
void LAPACK_clapmr(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    lapack_complex_float* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_dlapmr LAPACK_GLOBAL(dlapmr,DLAPMR)
void LAPACK_dlapmr(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    double* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_slapmr LAPACK_GLOBAL(slapmr,SLAPMR)
void LAPACK_slapmr(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    float* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_zlapmr LAPACK_GLOBAL(zlapmr,ZLAPMR)
void LAPACK_zlapmr(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    lapack_complex_double* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_clapmt LAPACK_GLOBAL(clapmt,CLAPMT)
void LAPACK_clapmt(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    lapack_complex_float* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_dlapmt LAPACK_GLOBAL(dlapmt,DLAPMT)
void LAPACK_dlapmt(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    double* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_slapmt LAPACK_GLOBAL(slapmt,SLAPMT)
void LAPACK_slapmt(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    float* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_zlapmt LAPACK_GLOBAL(zlapmt,ZLAPMT)
void LAPACK_zlapmt(
    lapack_logical const* forwrd, lapack_int const* m, lapack_int const* n,
    lapack_complex_double* X, lapack_int const* ldx, lapack_int* K );

#define LAPACK_dlapy2 LAPACK_GLOBAL(dlapy2,DLAPY2)
double LAPACK_dlapy2(
    double const* x,
    double const* y );

#define LAPACK_slapy2 LAPACK_GLOBAL(slapy2,SLAPY2)
lapack_float_return LAPACK_slapy2(
    float const* x,
    float const* y );

#define LAPACK_dlapy3 LAPACK_GLOBAL(dlapy3,DLAPY3)
double LAPACK_dlapy3(
    double const* x,
    double const* y,
    double const* z );

#define LAPACK_slapy3 LAPACK_GLOBAL(slapy3,SLAPY3)
lapack_float_return LAPACK_slapy3(
    float const* x,
    float const* y,
    float const* z );

#define LAPACK_clarcm LAPACK_GLOBAL(clarcm,CLARCM)
void LAPACK_clarcm(
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* C, lapack_int const* ldc,
    float* rwork );

#define LAPACK_zlarcm LAPACK_GLOBAL(zlarcm,ZLARCM)
void LAPACK_zlarcm(
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* C, lapack_int const* ldc,
    double* rwork );

#define LAPACK_clarf LAPACK_GLOBAL(clarf,CLARF)
void LAPACK_clarf(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* V, lapack_int const* incv,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work );

#define LAPACK_dlarf LAPACK_GLOBAL(dlarf,DLARF)
void LAPACK_dlarf(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    double const* V, lapack_int const* incv,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work );

#define LAPACK_slarf LAPACK_GLOBAL(slarf,SLARF)
void LAPACK_slarf(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    float const* V, lapack_int const* incv,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work );

#define LAPACK_zlarf LAPACK_GLOBAL(zlarf,ZLARF)
void LAPACK_zlarf(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* V, lapack_int const* incv,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work );

#define LAPACK_clarfb LAPACK_GLOBAL(clarfb,CLARFB)
void LAPACK_clarfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* V, lapack_int const* ldv,
    lapack_complex_float const* T, lapack_int const* ldt,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* ldwork );

#define LAPACK_dlarfb LAPACK_GLOBAL(dlarfb,DLARFB)
void LAPACK_dlarfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* V, lapack_int const* ldv,
    double const* T, lapack_int const* ldt,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* ldwork );

#define LAPACK_slarfb LAPACK_GLOBAL(slarfb,SLARFB)
void LAPACK_slarfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* V, lapack_int const* ldv,
    float const* T, lapack_int const* ldt,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* ldwork );

#define LAPACK_zlarfb LAPACK_GLOBAL(zlarfb,ZLARFB)
void LAPACK_zlarfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* V, lapack_int const* ldv,
    lapack_complex_double const* T, lapack_int const* ldt,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* ldwork );

#define LAPACK_clarfg LAPACK_GLOBAL(clarfg,CLARFG)
void LAPACK_clarfg(
    lapack_int const* n,
    lapack_complex_float* alpha,
    lapack_complex_float* X, lapack_int const* incx,
    lapack_complex_float* tau );

#define LAPACK_dlarfg LAPACK_GLOBAL(dlarfg,DLARFG)
void LAPACK_dlarfg(
    lapack_int const* n,
    double* alpha,
    double* X, lapack_int const* incx,
    double* tau );

#define LAPACK_slarfg LAPACK_GLOBAL(slarfg,SLARFG)
void LAPACK_slarfg(
    lapack_int const* n,
    float* alpha,
    float* X, lapack_int const* incx,
    float* tau );

#define LAPACK_zlarfg LAPACK_GLOBAL(zlarfg,ZLARFG)
void LAPACK_zlarfg(
    lapack_int const* n,
    lapack_complex_double* alpha,
    lapack_complex_double* X, lapack_int const* incx,
    lapack_complex_double* tau );

#define LAPACK_clarft LAPACK_GLOBAL(clarft,CLARFT)
void LAPACK_clarft(
    char const* direct, char const* storev,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* V, lapack_int const* ldv,
    lapack_complex_float const* tau,
    lapack_complex_float* T, lapack_int const* ldt );

#define LAPACK_dlarft LAPACK_GLOBAL(dlarft,DLARFT)
void LAPACK_dlarft(
    char const* direct, char const* storev,
    lapack_int const* n, lapack_int const* k,
    double const* V, lapack_int const* ldv,
    double const* tau,
    double* T, lapack_int const* ldt );

#define LAPACK_slarft LAPACK_GLOBAL(slarft,SLARFT)
void LAPACK_slarft(
    char const* direct, char const* storev,
    lapack_int const* n, lapack_int const* k,
    float const* V, lapack_int const* ldv,
    float const* tau,
    float* T, lapack_int const* ldt );

#define LAPACK_zlarft LAPACK_GLOBAL(zlarft,ZLARFT)
void LAPACK_zlarft(
    char const* direct, char const* storev,
    lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* V, lapack_int const* ldv,
    lapack_complex_double const* tau,
    lapack_complex_double* T, lapack_int const* ldt );

#define LAPACK_clarfx LAPACK_GLOBAL(clarfx,CLARFX)
void LAPACK_clarfx(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* V,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work );

#define LAPACK_dlarfx LAPACK_GLOBAL(dlarfx,DLARFX)
void LAPACK_dlarfx(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    double const* V,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work );

#define LAPACK_slarfx LAPACK_GLOBAL(slarfx,SLARFX)
void LAPACK_slarfx(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    float const* V,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work );

#define LAPACK_zlarfx LAPACK_GLOBAL(zlarfx,ZLARFX)
void LAPACK_zlarfx(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* V,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work );

#define LAPACK_clarnv LAPACK_GLOBAL(clarnv,CLARNV)
void LAPACK_clarnv(
    lapack_int const* idist, lapack_int* iseed, lapack_int const* n,
    lapack_complex_float* X );

#define LAPACK_dlarnv LAPACK_GLOBAL(dlarnv,DLARNV)
void LAPACK_dlarnv(
    lapack_int const* idist, lapack_int* iseed, lapack_int const* n,
    double* X );

#define LAPACK_slarnv LAPACK_GLOBAL(slarnv,SLARNV)
void LAPACK_slarnv(
    lapack_int const* idist, lapack_int* iseed, lapack_int const* n,
    float* X );

#define LAPACK_zlarnv LAPACK_GLOBAL(zlarnv,ZLARNV)
void LAPACK_zlarnv(
    lapack_int const* idist, lapack_int* iseed, lapack_int const* n,
    lapack_complex_double* X );

#define LAPACK_dlartgp LAPACK_GLOBAL(dlartgp,DLARTGP)
void LAPACK_dlartgp(
    double const* f,
    double const* g,
    double* cs,
    double* sn,
    double* r );

#define LAPACK_slartgp LAPACK_GLOBAL(slartgp,SLARTGP)
void LAPACK_slartgp(
    float const* f,
    float const* g,
    float* cs,
    float* sn,
    float* r );

#define LAPACK_dlartgs LAPACK_GLOBAL(dlartgs,DLARTGS)
void LAPACK_dlartgs(
    double const* x,
    double const* y,
    double const* sigma,
    double* cs,
    double* sn );

#define LAPACK_slartgs LAPACK_GLOBAL(slartgs,SLARTGS)
void LAPACK_slartgs(
    float const* x,
    float const* y,
    float const* sigma,
    float* cs,
    float* sn );

#define LAPACK_clascl LAPACK_GLOBAL(clascl,CLASCL)
void LAPACK_clascl(
    char const* type,
    lapack_int const* kl, lapack_int const* ku,
    float const* cfrom,
    float const* cto, lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dlascl LAPACK_GLOBAL(dlascl,DLASCL)
void LAPACK_dlascl(
    char const* type,
    lapack_int const* kl, lapack_int const* ku,
    double const* cfrom,
    double const* cto, lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_slascl LAPACK_GLOBAL(slascl,SLASCL)
void LAPACK_slascl(
    char const* type,
    lapack_int const* kl, lapack_int const* ku,
    float const* cfrom,
    float const* cto, lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_zlascl LAPACK_GLOBAL(zlascl,ZLASCL)
void LAPACK_zlascl(
    char const* type,
    lapack_int const* kl, lapack_int const* ku,
    double const* cfrom,
    double const* cto, lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_claset LAPACK_GLOBAL(claset,CLASET)
void LAPACK_claset(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* alpha,
    lapack_complex_float const* beta,
    lapack_complex_float* A, lapack_int const* lda );

#define LAPACK_dlaset LAPACK_GLOBAL(dlaset,DLASET)
void LAPACK_dlaset(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    double const* alpha,
    double const* beta,
    double* A, lapack_int const* lda );

#define LAPACK_slaset LAPACK_GLOBAL(slaset,SLASET)
void LAPACK_slaset(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    float const* alpha,
    float const* beta,
    float* A, lapack_int const* lda );

#define LAPACK_zlaset LAPACK_GLOBAL(zlaset,ZLASET)
void LAPACK_zlaset(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* alpha,
    lapack_complex_double const* beta,
    lapack_complex_double* A, lapack_int const* lda );

#define LAPACK_dlasrt LAPACK_GLOBAL(dlasrt,DLASRT)
void LAPACK_dlasrt(
    char const* id,
    lapack_int const* n,
    double* D,
    lapack_int* info );

#define LAPACK_slasrt LAPACK_GLOBAL(slasrt,SLASRT)
void LAPACK_slasrt(
    char const* id,
    lapack_int const* n,
    float* D,
    lapack_int* info );

#define LAPACK_classq LAPACK_GLOBAL(classq,CLASSQ)
void LAPACK_classq(
    lapack_int const* n,
    lapack_complex_float const* X, lapack_int const* incx,
    float* scale,
    float* sumsq );

#define LAPACK_dlassq LAPACK_GLOBAL(dlassq,DLASSQ)
void LAPACK_dlassq(
    lapack_int const* n,
    double const* X, lapack_int const* incx,
    double* scale,
    double* sumsq );

#define LAPACK_slassq LAPACK_GLOBAL(slassq,SLASSQ)
void LAPACK_slassq(
    lapack_int const* n,
    float const* X, lapack_int const* incx,
    float* scale,
    float* sumsq );

#define LAPACK_zlassq LAPACK_GLOBAL(zlassq,ZLASSQ)
void LAPACK_zlassq(
    lapack_int const* n,
    lapack_complex_double const* X, lapack_int const* incx,
    double* scale,
    double* sumsq );

#define LAPACK_claswp LAPACK_GLOBAL(claswp,CLASWP)
void LAPACK_claswp(
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* k1, lapack_int const* k2, lapack_int const* ipiv, lapack_int const* incx );

#define LAPACK_dlaswp LAPACK_GLOBAL(dlaswp,DLASWP)
void LAPACK_dlaswp(
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int const* k1, lapack_int const* k2, lapack_int const* ipiv, lapack_int const* incx );

#define LAPACK_slaswp LAPACK_GLOBAL(slaswp,SLASWP)
void LAPACK_slaswp(
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int const* k1, lapack_int const* k2, lapack_int const* ipiv, lapack_int const* incx );

#define LAPACK_zlaswp LAPACK_GLOBAL(zlaswp,ZLASWP)
void LAPACK_zlaswp(
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* k1, lapack_int const* k2, lapack_int const* ipiv, lapack_int const* incx );

#define LAPACK_clatms LAPACK_GLOBAL(clatms,CLATMS)
void LAPACK_clatms(
    lapack_int const* m, lapack_int const* n, char const* dist,
    lapack_int* iseed, char const* sym,
    float* D,
    lapack_int const* mode,
    float const* cond,
    float const* dmax, lapack_int const* kl, lapack_int const* ku, char const* pack,
    lapack_complex_float* A,
    lapack_int const* lda,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dlatms LAPACK_GLOBAL(dlatms,DLATMS)
void LAPACK_dlatms(
    lapack_int const* m, lapack_int const* n, char const* dist,
    lapack_int* iseed, char const* sym,
    double* D,
    lapack_int const* mode,
    double const* cond,
    double const* dmax, lapack_int const* kl, lapack_int const* ku, char const* pack,
    double* A,
    lapack_int const* lda,
    double* work,
    lapack_int* info );

#define LAPACK_slatms LAPACK_GLOBAL(slatms,SLATMS)
void LAPACK_slatms(
    lapack_int const* m, lapack_int const* n, char const* dist,
    lapack_int* iseed, char const* sym,
    float* D,
    lapack_int const* mode,
    float const* cond,
    float const* dmax, lapack_int const* kl, lapack_int const* ku, char const* pack,
    float* A,
    lapack_int const* lda,
    float* work,
    lapack_int* info );

#define LAPACK_zlatms LAPACK_GLOBAL(zlatms,ZLATMS)
void LAPACK_zlatms(
    lapack_int const* m, lapack_int const* n, char const* dist,
    lapack_int* iseed, char const* sym,
    double* D,
    lapack_int const* mode,
    double const* cond,
    double const* dmax, lapack_int const* kl, lapack_int const* ku, char const* pack,
    lapack_complex_double* A,
    lapack_int const* lda,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_clauum LAPACK_GLOBAL(clauum,CLAUUM)
void LAPACK_clauum(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dlauum LAPACK_GLOBAL(dlauum,DLAUUM)
void LAPACK_dlauum(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_slauum LAPACK_GLOBAL(slauum,SLAUUM)
void LAPACK_slauum(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_zlauum LAPACK_GLOBAL(zlauum,ZLAUUM)
void LAPACK_zlauum(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_ilaver LAPACK_GLOBAL(ilaver,ILAVER)
void LAPACK_ilaver(
    lapack_int* vers_major, lapack_int* vers_minor, lapack_int* vers_patch );

#define LAPACK_dopgtr LAPACK_GLOBAL(dopgtr,DOPGTR)
void LAPACK_dopgtr(
    char const* uplo,
    lapack_int const* n,
    double const* AP,
    double const* tau,
    double* Q, lapack_int const* ldq,
    double* work,
    lapack_int* info );

#define LAPACK_sopgtr LAPACK_GLOBAL(sopgtr,SOPGTR)
void LAPACK_sopgtr(
    char const* uplo,
    lapack_int const* n,
    float const* AP,
    float const* tau,
    float* Q, lapack_int const* ldq,
    float* work,
    lapack_int* info );

#define LAPACK_dopmtr LAPACK_GLOBAL(dopmtr,DOPMTR)
void LAPACK_dopmtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    double const* AP,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work,
    lapack_int* info );

#define LAPACK_sopmtr LAPACK_GLOBAL(sopmtr,SOPMTR)
void LAPACK_sopmtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    float const* AP,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work,
    lapack_int* info );

#define LAPACK_dorbdb LAPACK_GLOBAL(dorbdb,DORBDB)
void LAPACK_dorbdb(
    char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    double* X11, lapack_int const* ldx11,
    double* X12, lapack_int const* ldx12,
    double* X21, lapack_int const* ldx21,
    double* X22, lapack_int const* ldx22,
    double* theta,
    double* phi,
    double* TAUP1,
    double* TAUP2,
    double* TAUQ1,
    double* TAUQ2,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorbdb LAPACK_GLOBAL(sorbdb,SORBDB)
void LAPACK_sorbdb(
    char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    float* X11, lapack_int const* ldx11,
    float* X12, lapack_int const* ldx12,
    float* X21, lapack_int const* ldx21,
    float* X22, lapack_int const* ldx22,
    float* theta,
    float* phi,
    float* TAUP1,
    float* TAUP2,
    float* TAUQ1,
    float* TAUQ2,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dorcsd LAPACK_GLOBAL(dorcsd,DORCSD)
void LAPACK_dorcsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    double* X11, lapack_int const* ldx11,
    double* X12, lapack_int const* ldx12,
    double* X21, lapack_int const* ldx21,
    double* X22, lapack_int const* ldx22,
    double* theta,
    double* U1, lapack_int const* ldu1,
    double* U2, lapack_int const* ldu2,
    double* V1T, lapack_int const* ldv1t,
    double* V2T, lapack_int const* ldv2t,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sorcsd LAPACK_GLOBAL(sorcsd,SORCSD)
void LAPACK_sorcsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    float* X11, lapack_int const* ldx11,
    float* X12, lapack_int const* ldx12,
    float* X21, lapack_int const* ldx21,
    float* X22, lapack_int const* ldx22,
    float* theta,
    float* U1, lapack_int const* ldu1,
    float* U2, lapack_int const* ldu2,
    float* V1T, lapack_int const* ldv1t,
    float* V2T, lapack_int const* ldv2t,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dorcsd2by1 LAPACK_GLOBAL(dorcsd2by1,DORCSD2BY1)
void LAPACK_dorcsd2by1(
    char const* jobu1, char const* jobu2, char const* jobv1t,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    double* X11, lapack_int const* ldx11,
    double* X21, lapack_int const* ldx21,
    double* theta,
    double* U1, lapack_int const* ldu1,
    double* U2, lapack_int const* ldu2,
    double* V1T, lapack_int const* ldv1t,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sorcsd2by1 LAPACK_GLOBAL(sorcsd2by1,SORCSD2BY1)
void LAPACK_sorcsd2by1(
    char const* jobu1, char const* jobu2, char const* jobv1t,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    float* X11, lapack_int const* ldx11,
    float* X21, lapack_int const* ldx21,
    float* theta,
    float* U1, lapack_int const* ldu1,
    float* U2, lapack_int const* ldu2,
    float* V1T, lapack_int const* ldv1t,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dorgbr LAPACK_GLOBAL(dorgbr,DORGBR)
void LAPACK_dorgbr(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorgbr LAPACK_GLOBAL(sorgbr,SORGBR)
void LAPACK_sorgbr(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float* A, lapack_int const* lda,
    float const* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dorghr LAPACK_GLOBAL(dorghr,DORGHR)
void LAPACK_dorghr(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorghr LAPACK_GLOBAL(sorghr,SORGHR)
void LAPACK_sorghr(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float* A, lapack_int const* lda,
    float const* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dorglq LAPACK_GLOBAL(dorglq,DORGLQ)
void LAPACK_dorglq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorglq LAPACK_GLOBAL(sorglq,SORGLQ)
void LAPACK_sorglq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float* A, lapack_int const* lda,
    float const* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dorgql LAPACK_GLOBAL(dorgql,DORGQL)
void LAPACK_dorgql(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorgql LAPACK_GLOBAL(sorgql,SORGQL)
void LAPACK_sorgql(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float* A, lapack_int const* lda,
    float const* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dorgqr LAPACK_GLOBAL(dorgqr,DORGQR)
void LAPACK_dorgqr(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorgqr LAPACK_GLOBAL(sorgqr,SORGQR)
void LAPACK_sorgqr(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float* A, lapack_int const* lda,
    float const* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dorgrq LAPACK_GLOBAL(dorgrq,DORGRQ)
void LAPACK_dorgrq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorgrq LAPACK_GLOBAL(sorgrq,SORGRQ)
void LAPACK_sorgrq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float* A, lapack_int const* lda,
    float const* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dorgtr LAPACK_GLOBAL(dorgtr,DORGTR)
void LAPACK_dorgtr(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sorgtr LAPACK_GLOBAL(sorgtr,SORGTR)
void LAPACK_sorgtr(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormbr LAPACK_GLOBAL(dormbr,DORMBR)
void LAPACK_dormbr(
    char const* vect, char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormbr LAPACK_GLOBAL(sormbr,SORMBR)
void LAPACK_sormbr(
    char const* vect, char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormhr LAPACK_GLOBAL(dormhr,DORMHR)
void LAPACK_dormhr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormhr LAPACK_GLOBAL(sormhr,SORMHR)
void LAPACK_sormhr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormlq LAPACK_GLOBAL(dormlq,DORMLQ)
void LAPACK_dormlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormlq LAPACK_GLOBAL(sormlq,SORMLQ)
void LAPACK_sormlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormql LAPACK_GLOBAL(dormql,DORMQL)
void LAPACK_dormql(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormql LAPACK_GLOBAL(sormql,SORMQL)
void LAPACK_sormql(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormqr LAPACK_GLOBAL(dormqr,DORMQR)
void LAPACK_dormqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormqr LAPACK_GLOBAL(sormqr,SORMQR)
void LAPACK_sormqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormrq LAPACK_GLOBAL(dormrq,DORMRQ)
void LAPACK_dormrq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormrq LAPACK_GLOBAL(sormrq,SORMRQ)
void LAPACK_sormrq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormrz LAPACK_GLOBAL(dormrz,DORMRZ)
void LAPACK_dormrz(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormrz LAPACK_GLOBAL(sormrz,SORMRZ)
void LAPACK_sormrz(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dormtr LAPACK_GLOBAL(dormtr,DORMTR)
void LAPACK_dormtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_sormtr LAPACK_GLOBAL(sormtr,SORMTR)
void LAPACK_sormtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float const* tau,
    float* C, lapack_int const* ldc,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cpbcon LAPACK_GLOBAL(cpbcon,CPBCON)
void LAPACK_cpbcon(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dpbcon LAPACK_GLOBAL(dpbcon,DPBCON)
void LAPACK_dpbcon(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double const* AB, lapack_int const* ldab,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_spbcon LAPACK_GLOBAL(spbcon,SPBCON)
void LAPACK_spbcon(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float const* AB, lapack_int const* ldab,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zpbcon LAPACK_GLOBAL(zpbcon,ZPBCON)
void LAPACK_zpbcon(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpbequ LAPACK_GLOBAL(cpbequ,CPBEQU)
void LAPACK_cpbequ(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_dpbequ LAPACK_GLOBAL(dpbequ,DPBEQU)
void LAPACK_dpbequ(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double const* AB, lapack_int const* ldab,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_spbequ LAPACK_GLOBAL(spbequ,SPBEQU)
void LAPACK_spbequ(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float const* AB, lapack_int const* ldab,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_zpbequ LAPACK_GLOBAL(zpbequ,ZPBEQU)
void LAPACK_zpbequ(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_cpbrfs LAPACK_GLOBAL(cpbrfs,CPBRFS)
void LAPACK_cpbrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_float const* AB, lapack_int const* ldab,
    lapack_complex_float const* AFB, lapack_int const* ldafb,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dpbrfs LAPACK_GLOBAL(dpbrfs,DPBRFS)
void LAPACK_dpbrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab,
    double const* AFB, lapack_int const* ldafb,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_spbrfs LAPACK_GLOBAL(spbrfs,SPBRFS)
void LAPACK_spbrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    float const* AB, lapack_int const* ldab,
    float const* AFB, lapack_int const* ldafb,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zpbrfs LAPACK_GLOBAL(zpbrfs,ZPBRFS)
void LAPACK_zpbrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_double const* AB, lapack_int const* ldab,
    lapack_complex_double const* AFB, lapack_int const* ldafb,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpbstf LAPACK_GLOBAL(cpbstf,CPBSTF)
void LAPACK_cpbstf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_dpbstf LAPACK_GLOBAL(dpbstf,DPBSTF)
void LAPACK_dpbstf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_spbstf LAPACK_GLOBAL(spbstf,SPBSTF)
void LAPACK_spbstf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_zpbstf LAPACK_GLOBAL(zpbstf,ZPBSTF)
void LAPACK_zpbstf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_cpbsv LAPACK_GLOBAL(cpbsv,CPBSV)
void LAPACK_cpbsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dpbsv LAPACK_GLOBAL(dpbsv,DPBSV)
void LAPACK_dpbsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    double* AB, lapack_int const* ldab,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_spbsv LAPACK_GLOBAL(spbsv,SPBSV)
void LAPACK_spbsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    float* AB, lapack_int const* ldab,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zpbsv LAPACK_GLOBAL(zpbsv,ZPBSV)
void LAPACK_zpbsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cpbsvx LAPACK_GLOBAL(cpbsvx,CPBSVX)
void LAPACK_cpbsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_complex_float* AFB, lapack_int const* ldafb, char* equed,
    float* S,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dpbsvx LAPACK_GLOBAL(dpbsvx,DPBSVX)
void LAPACK_dpbsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    double* AB, lapack_int const* ldab,
    double* AFB, lapack_int const* ldafb, char* equed,
    double* S,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_spbsvx LAPACK_GLOBAL(spbsvx,SPBSVX)
void LAPACK_spbsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    float* AB, lapack_int const* ldab,
    float* AFB, lapack_int const* ldafb, char* equed,
    float* S,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zpbsvx LAPACK_GLOBAL(zpbsvx,ZPBSVX)
void LAPACK_zpbsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_complex_double* AFB, lapack_int const* ldafb, char* equed,
    double* S,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpbtrf LAPACK_GLOBAL(cpbtrf,CPBTRF)
void LAPACK_cpbtrf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_dpbtrf LAPACK_GLOBAL(dpbtrf,DPBTRF)
void LAPACK_dpbtrf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_spbtrf LAPACK_GLOBAL(spbtrf,SPBTRF)
void LAPACK_spbtrf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_zpbtrf LAPACK_GLOBAL(zpbtrf,ZPBTRF)
void LAPACK_zpbtrf(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double* AB, lapack_int const* ldab,
    lapack_int* info );

#define LAPACK_cpbtrs LAPACK_GLOBAL(cpbtrs,CPBTRS)
void LAPACK_cpbtrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_float const* AB, lapack_int const* ldab,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dpbtrs LAPACK_GLOBAL(dpbtrs,DPBTRS)
void LAPACK_dpbtrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_spbtrs LAPACK_GLOBAL(spbtrs,SPBTRS)
void LAPACK_spbtrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    float const* AB, lapack_int const* ldab,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zpbtrs LAPACK_GLOBAL(zpbtrs,ZPBTRS)
void LAPACK_zpbtrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_double const* AB, lapack_int const* ldab,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cpftrf LAPACK_GLOBAL(cpftrf,CPFTRF)
void LAPACK_cpftrf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A,
    lapack_int* info );

#define LAPACK_dpftrf LAPACK_GLOBAL(dpftrf,DPFTRF)
void LAPACK_dpftrf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    double* A,
    lapack_int* info );

#define LAPACK_spftrf LAPACK_GLOBAL(spftrf,SPFTRF)
void LAPACK_spftrf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    float* A,
    lapack_int* info );

#define LAPACK_zpftrf LAPACK_GLOBAL(zpftrf,ZPFTRF)
void LAPACK_zpftrf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A,
    lapack_int* info );

#define LAPACK_cpftri LAPACK_GLOBAL(cpftri,CPFTRI)
void LAPACK_cpftri(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A,
    lapack_int* info );

#define LAPACK_dpftri LAPACK_GLOBAL(dpftri,DPFTRI)
void LAPACK_dpftri(
    char const* transr, char const* uplo,
    lapack_int const* n,
    double* A,
    lapack_int* info );

#define LAPACK_spftri LAPACK_GLOBAL(spftri,SPFTRI)
void LAPACK_spftri(
    char const* transr, char const* uplo,
    lapack_int const* n,
    float* A,
    lapack_int* info );

#define LAPACK_zpftri LAPACK_GLOBAL(zpftri,ZPFTRI)
void LAPACK_zpftri(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A,
    lapack_int* info );

#define LAPACK_cpftrs LAPACK_GLOBAL(cpftrs,CPFTRS)
void LAPACK_cpftrs(
    char const* transr, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dpftrs LAPACK_GLOBAL(dpftrs,DPFTRS)
void LAPACK_dpftrs(
    char const* transr, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_spftrs LAPACK_GLOBAL(spftrs,SPFTRS)
void LAPACK_spftrs(
    char const* transr, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zpftrs LAPACK_GLOBAL(zpftrs,ZPFTRS)
void LAPACK_zpftrs(
    char const* transr, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cpocon LAPACK_GLOBAL(cpocon,CPOCON)
void LAPACK_cpocon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dpocon LAPACK_GLOBAL(dpocon,DPOCON)
void LAPACK_dpocon(
    char const* uplo,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_spocon LAPACK_GLOBAL(spocon,SPOCON)
void LAPACK_spocon(
    char const* uplo,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zpocon LAPACK_GLOBAL(zpocon,ZPOCON)
void LAPACK_zpocon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpoequ LAPACK_GLOBAL(cpoequ,CPOEQU)
void LAPACK_cpoequ(
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_dpoequ LAPACK_GLOBAL(dpoequ,DPOEQU)
void LAPACK_dpoequ(
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_spoequ LAPACK_GLOBAL(spoequ,SPOEQU)
void LAPACK_spoequ(
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_zpoequ LAPACK_GLOBAL(zpoequ,ZPOEQU)
void LAPACK_zpoequ(
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_cpoequb LAPACK_GLOBAL(cpoequb,CPOEQUB)
void LAPACK_cpoequb(
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_dpoequb LAPACK_GLOBAL(dpoequb,DPOEQUB)
void LAPACK_dpoequb(
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_spoequb LAPACK_GLOBAL(spoequb,SPOEQUB)
void LAPACK_spoequb(
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_zpoequb LAPACK_GLOBAL(zpoequb,ZPOEQUB)
void LAPACK_zpoequb(
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_cporfs LAPACK_GLOBAL(cporfs,CPORFS)
void LAPACK_cporfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dporfs LAPACK_GLOBAL(dporfs,DPORFS)
void LAPACK_dporfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* AF, lapack_int const* ldaf,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sporfs LAPACK_GLOBAL(sporfs,SPORFS)
void LAPACK_sporfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* AF, lapack_int const* ldaf,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zporfs LAPACK_GLOBAL(zporfs,ZPORFS)
void LAPACK_zporfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cporfsx LAPACK_GLOBAL(cporfsx,CPORFSX)
void LAPACK_cporfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf,
    float* S,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dporfsx LAPACK_GLOBAL(dporfsx,DPORFSX)
void LAPACK_dporfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* AF, lapack_int const* ldaf,
    double* S,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sporfsx LAPACK_GLOBAL(sporfsx,SPORFSX)
void LAPACK_sporfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* AF, lapack_int const* ldaf,
    float* S,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zporfsx LAPACK_GLOBAL(zporfsx,ZPORFSX)
void LAPACK_zporfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf,
    double* S,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cposv LAPACK_GLOBAL(cposv,CPOSV)
void LAPACK_cposv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dposv LAPACK_GLOBAL(dposv,DPOSV)
void LAPACK_dposv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sposv LAPACK_GLOBAL(sposv,SPOSV)
void LAPACK_sposv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zposv LAPACK_GLOBAL(zposv,ZPOSV)
void LAPACK_zposv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsposv LAPACK_GLOBAL(dsposv,DSPOSV)
void LAPACK_dsposv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* work,
    float* swork, lapack_int* iter,
    lapack_int* info );

#define LAPACK_zcposv LAPACK_GLOBAL(zcposv,ZCPOSV)
void LAPACK_zcposv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    lapack_complex_double* work,
    lapack_complex_float* swork,
    double* rwork, lapack_int* iter,
    lapack_int* info );

#define LAPACK_cposvx LAPACK_GLOBAL(cposvx,CPOSVX)
void LAPACK_cposvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, char* equed,
    float* S,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dposvx LAPACK_GLOBAL(dposvx,DPOSVX)
void LAPACK_dposvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* AF, lapack_int const* ldaf, char* equed,
    double* S,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sposvx LAPACK_GLOBAL(sposvx,SPOSVX)
void LAPACK_sposvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* AF, lapack_int const* ldaf, char* equed,
    float* S,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zposvx LAPACK_GLOBAL(zposvx,ZPOSVX)
void LAPACK_zposvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, char* equed,
    double* S,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cposvxx LAPACK_GLOBAL(cposvxx,CPOSVXX)
void LAPACK_cposvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, char* equed,
    float* S,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dposvxx LAPACK_GLOBAL(dposvxx,DPOSVXX)
void LAPACK_dposvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* AF, lapack_int const* ldaf, char* equed,
    double* S,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sposvxx LAPACK_GLOBAL(sposvxx,SPOSVXX)
void LAPACK_sposvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* AF, lapack_int const* ldaf, char* equed,
    float* S,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zposvxx LAPACK_GLOBAL(zposvxx,ZPOSVXX)
void LAPACK_zposvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, char* equed,
    double* S,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpotf2 LAPACK_GLOBAL(cpotf2,CPOTF2)
void LAPACK_cpotf2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dpotf2 LAPACK_GLOBAL(dpotf2,DPOTF2)
void LAPACK_dpotf2(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_spotf2 LAPACK_GLOBAL(spotf2,SPOTF2)
void LAPACK_spotf2(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_zpotf2 LAPACK_GLOBAL(zpotf2,ZPOTF2)
void LAPACK_zpotf2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_cpotrf LAPACK_GLOBAL(cpotrf,CPOTRF)
void LAPACK_cpotrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dpotrf LAPACK_GLOBAL(dpotrf,DPOTRF)
void LAPACK_dpotrf(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_spotrf LAPACK_GLOBAL(spotrf,SPOTRF)
void LAPACK_spotrf(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_zpotrf LAPACK_GLOBAL(zpotrf,ZPOTRF)
void LAPACK_zpotrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_cpotrf2 LAPACK_GLOBAL(cpotrf2,CPOTRF2)
void LAPACK_cpotrf2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dpotrf2 LAPACK_GLOBAL(dpotrf2,DPOTRF2)
void LAPACK_dpotrf2(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_spotrf2 LAPACK_GLOBAL(spotrf2,SPOTRF2)
void LAPACK_spotrf2(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_zpotrf2 LAPACK_GLOBAL(zpotrf2,ZPOTRF2)
void LAPACK_zpotrf2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_cpotri LAPACK_GLOBAL(cpotri,CPOTRI)
void LAPACK_cpotri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dpotri LAPACK_GLOBAL(dpotri,DPOTRI)
void LAPACK_dpotri(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_spotri LAPACK_GLOBAL(spotri,SPOTRI)
void LAPACK_spotri(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_zpotri LAPACK_GLOBAL(zpotri,ZPOTRI)
void LAPACK_zpotri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_cpotrs LAPACK_GLOBAL(cpotrs,CPOTRS)
void LAPACK_cpotrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dpotrs LAPACK_GLOBAL(dpotrs,DPOTRS)
void LAPACK_dpotrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_spotrs LAPACK_GLOBAL(spotrs,SPOTRS)
void LAPACK_spotrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zpotrs LAPACK_GLOBAL(zpotrs,ZPOTRS)
void LAPACK_zpotrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cppcon LAPACK_GLOBAL(cppcon,CPPCON)
void LAPACK_cppcon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dppcon LAPACK_GLOBAL(dppcon,DPPCON)
void LAPACK_dppcon(
    char const* uplo,
    lapack_int const* n,
    double const* AP,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sppcon LAPACK_GLOBAL(sppcon,SPPCON)
void LAPACK_sppcon(
    char const* uplo,
    lapack_int const* n,
    float const* AP,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zppcon LAPACK_GLOBAL(zppcon,ZPPCON)
void LAPACK_zppcon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cppequ LAPACK_GLOBAL(cppequ,CPPEQU)
void LAPACK_cppequ(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_dppequ LAPACK_GLOBAL(dppequ,DPPEQU)
void LAPACK_dppequ(
    char const* uplo,
    lapack_int const* n,
    double const* AP,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_sppequ LAPACK_GLOBAL(sppequ,SPPEQU)
void LAPACK_sppequ(
    char const* uplo,
    lapack_int const* n,
    float const* AP,
    float* S,
    float* scond,
    float* amax,
    lapack_int* info );

#define LAPACK_zppequ LAPACK_GLOBAL(zppequ,ZPPEQU)
void LAPACK_zppequ(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP,
    double* S,
    double* scond,
    double* amax,
    lapack_int* info );

#define LAPACK_cpprfs LAPACK_GLOBAL(cpprfs,CPPRFS)
void LAPACK_cpprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float const* AFP,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dpprfs LAPACK_GLOBAL(dpprfs,DPPRFS)
void LAPACK_dpprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* AP,
    double const* AFP,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_spprfs LAPACK_GLOBAL(spprfs,SPPRFS)
void LAPACK_spprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* AP,
    float const* AFP,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zpprfs LAPACK_GLOBAL(zpprfs,ZPPRFS)
void LAPACK_zpprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double const* AFP,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cppsv LAPACK_GLOBAL(cppsv,CPPSV)
void LAPACK_cppsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* AP,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dppsv LAPACK_GLOBAL(dppsv,DPPSV)
void LAPACK_dppsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* AP,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sppsv LAPACK_GLOBAL(sppsv,SPPSV)
void LAPACK_sppsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* AP,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zppsv LAPACK_GLOBAL(zppsv,ZPPSV)
void LAPACK_zppsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* AP,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cppsvx LAPACK_GLOBAL(cppsvx,CPPSVX)
void LAPACK_cppsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* AP,
    lapack_complex_float* AFP, char* equed,
    float* S,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dppsvx LAPACK_GLOBAL(dppsvx,DPPSVX)
void LAPACK_dppsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* AP,
    double* AFP, char* equed,
    double* S,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sppsvx LAPACK_GLOBAL(sppsvx,SPPSVX)
void LAPACK_sppsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* AP,
    float* AFP, char* equed,
    float* S,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zppsvx LAPACK_GLOBAL(zppsvx,ZPPSVX)
void LAPACK_zppsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* AP,
    lapack_complex_double* AFP, char* equed,
    double* S,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpptrf LAPACK_GLOBAL(cpptrf,CPPTRF)
void LAPACK_cpptrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    lapack_int* info );

#define LAPACK_dpptrf LAPACK_GLOBAL(dpptrf,DPPTRF)
void LAPACK_dpptrf(
    char const* uplo,
    lapack_int const* n,
    double* AP,
    lapack_int* info );

#define LAPACK_spptrf LAPACK_GLOBAL(spptrf,SPPTRF)
void LAPACK_spptrf(
    char const* uplo,
    lapack_int const* n,
    float* AP,
    lapack_int* info );

#define LAPACK_zpptrf LAPACK_GLOBAL(zpptrf,ZPPTRF)
void LAPACK_zpptrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    lapack_int* info );

#define LAPACK_cpptri LAPACK_GLOBAL(cpptri,CPPTRI)
void LAPACK_cpptri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP,
    lapack_int* info );

#define LAPACK_dpptri LAPACK_GLOBAL(dpptri,DPPTRI)
void LAPACK_dpptri(
    char const* uplo,
    lapack_int const* n,
    double* AP,
    lapack_int* info );

#define LAPACK_spptri LAPACK_GLOBAL(spptri,SPPTRI)
void LAPACK_spptri(
    char const* uplo,
    lapack_int const* n,
    float* AP,
    lapack_int* info );

#define LAPACK_zpptri LAPACK_GLOBAL(zpptri,ZPPTRI)
void LAPACK_zpptri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    lapack_int* info );

#define LAPACK_cpptrs LAPACK_GLOBAL(cpptrs,CPPTRS)
void LAPACK_cpptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dpptrs LAPACK_GLOBAL(dpptrs,DPPTRS)
void LAPACK_dpptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* AP,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_spptrs LAPACK_GLOBAL(spptrs,SPPTRS)
void LAPACK_spptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* AP,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zpptrs LAPACK_GLOBAL(zpptrs,ZPPTRS)
void LAPACK_zpptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cpstrf LAPACK_GLOBAL(cpstrf,CPSTRF)
void LAPACK_cpstrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* piv, lapack_int* rank,
    float const* tol,
    float* work,
    lapack_int* info );

#define LAPACK_dpstrf LAPACK_GLOBAL(dpstrf,DPSTRF)
void LAPACK_dpstrf(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* piv, lapack_int* rank,
    double const* tol,
    double* work,
    lapack_int* info );

#define LAPACK_spstrf LAPACK_GLOBAL(spstrf,SPSTRF)
void LAPACK_spstrf(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* piv, lapack_int* rank,
    float const* tol,
    float* work,
    lapack_int* info );

#define LAPACK_zpstrf LAPACK_GLOBAL(zpstrf,ZPSTRF)
void LAPACK_zpstrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* piv, lapack_int* rank,
    double const* tol,
    double* work,
    lapack_int* info );

#define LAPACK_cptcon LAPACK_GLOBAL(cptcon,CPTCON)
void LAPACK_cptcon(
    lapack_int const* n,
    float const* D,
    lapack_complex_float const* E,
    float const* anorm,
    float* rcond,
    float* rwork,
    lapack_int* info );

#define LAPACK_dptcon LAPACK_GLOBAL(dptcon,DPTCON)
void LAPACK_dptcon(
    lapack_int const* n,
    double const* D,
    double const* E,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* info );

#define LAPACK_sptcon LAPACK_GLOBAL(sptcon,SPTCON)
void LAPACK_sptcon(
    lapack_int const* n,
    float const* D,
    float const* E,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* info );

#define LAPACK_zptcon LAPACK_GLOBAL(zptcon,ZPTCON)
void LAPACK_zptcon(
    lapack_int const* n,
    double const* D,
    lapack_complex_double const* E,
    double const* anorm,
    double* rcond,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpteqr LAPACK_GLOBAL(cpteqr,CPTEQR)
void LAPACK_cpteqr(
    char const* compz,
    lapack_int const* n,
    float* D,
    float* E,
    lapack_complex_float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_dpteqr LAPACK_GLOBAL(dpteqr,DPTEQR)
void LAPACK_dpteqr(
    char const* compz,
    lapack_int const* n,
    double* D,
    double* E,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_spteqr LAPACK_GLOBAL(spteqr,SPTEQR)
void LAPACK_spteqr(
    char const* compz,
    lapack_int const* n,
    float* D,
    float* E,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_zpteqr LAPACK_GLOBAL(zpteqr,ZPTEQR)
void LAPACK_zpteqr(
    char const* compz,
    lapack_int const* n,
    double* D,
    double* E,
    lapack_complex_double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_cptrfs LAPACK_GLOBAL(cptrfs,CPTRFS)
void LAPACK_cptrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* D,
    lapack_complex_float const* E,
    float const* DF,
    lapack_complex_float const* EF,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dptrfs LAPACK_GLOBAL(dptrfs,DPTRFS)
void LAPACK_dptrfs(
    lapack_int const* n, lapack_int const* nrhs,
    double const* D,
    double const* E,
    double const* DF,
    double const* EF,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* info );

#define LAPACK_sptrfs LAPACK_GLOBAL(sptrfs,SPTRFS)
void LAPACK_sptrfs(
    lapack_int const* n, lapack_int const* nrhs,
    float const* D,
    float const* E,
    float const* DF,
    float const* EF,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* info );

#define LAPACK_zptrfs LAPACK_GLOBAL(zptrfs,ZPTRFS)
void LAPACK_zptrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* D,
    lapack_complex_double const* E,
    double const* DF,
    lapack_complex_double const* EF,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cptsv LAPACK_GLOBAL(cptsv,CPTSV)
void LAPACK_cptsv(
    lapack_int const* n, lapack_int const* nrhs,
    float* D,
    lapack_complex_float* E,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dptsv LAPACK_GLOBAL(dptsv,DPTSV)
void LAPACK_dptsv(
    lapack_int const* n, lapack_int const* nrhs,
    double* D,
    double* E,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sptsv LAPACK_GLOBAL(sptsv,SPTSV)
void LAPACK_sptsv(
    lapack_int const* n, lapack_int const* nrhs,
    float* D,
    float* E,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zptsv LAPACK_GLOBAL(zptsv,ZPTSV)
void LAPACK_zptsv(
    lapack_int const* n, lapack_int const* nrhs,
    double* D,
    lapack_complex_double* E,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cptsvx LAPACK_GLOBAL(cptsvx,CPTSVX)
void LAPACK_cptsvx(
    char const* fact,
    lapack_int const* n, lapack_int const* nrhs,
    float const* D,
    lapack_complex_float const* E,
    float* DF,
    lapack_complex_float* EF,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dptsvx LAPACK_GLOBAL(dptsvx,DPTSVX)
void LAPACK_dptsvx(
    char const* fact,
    lapack_int const* n, lapack_int const* nrhs,
    double const* D,
    double const* E,
    double* DF,
    double* EF,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* info );

#define LAPACK_sptsvx LAPACK_GLOBAL(sptsvx,SPTSVX)
void LAPACK_sptsvx(
    char const* fact,
    lapack_int const* n, lapack_int const* nrhs,
    float const* D,
    float const* E,
    float* DF,
    float* EF,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* info );

#define LAPACK_zptsvx LAPACK_GLOBAL(zptsvx,ZPTSVX)
void LAPACK_zptsvx(
    char const* fact,
    lapack_int const* n, lapack_int const* nrhs,
    double const* D,
    lapack_complex_double const* E,
    double* DF,
    lapack_complex_double* EF,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cpttrf LAPACK_GLOBAL(cpttrf,CPTTRF)
void LAPACK_cpttrf(
    lapack_int const* n,
    float* D,
    lapack_complex_float* E,
    lapack_int* info );

#define LAPACK_dpttrf LAPACK_GLOBAL(dpttrf,DPTTRF)
void LAPACK_dpttrf(
    lapack_int const* n,
    double* D,
    double* E,
    lapack_int* info );

#define LAPACK_spttrf LAPACK_GLOBAL(spttrf,SPTTRF)
void LAPACK_spttrf(
    lapack_int const* n,
    float* D,
    float* E,
    lapack_int* info );

#define LAPACK_zpttrf LAPACK_GLOBAL(zpttrf,ZPTTRF)
void LAPACK_zpttrf(
    lapack_int const* n,
    double* D,
    lapack_complex_double* E,
    lapack_int* info );

#define LAPACK_cpttrs LAPACK_GLOBAL(cpttrs,CPTTRS)
void LAPACK_cpttrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* D,
    lapack_complex_float const* E,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dpttrs LAPACK_GLOBAL(dpttrs,DPTTRS)
void LAPACK_dpttrs(
    lapack_int const* n, lapack_int const* nrhs,
    double const* D,
    double const* E,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_spttrs LAPACK_GLOBAL(spttrs,SPTTRS)
void LAPACK_spttrs(
    lapack_int const* n, lapack_int const* nrhs,
    float const* D,
    float const* E,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zpttrs LAPACK_GLOBAL(zpttrs,ZPTTRS)
void LAPACK_zpttrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* D,
    lapack_complex_double const* E,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsbev LAPACK_GLOBAL(dsbev,DSBEV)
void LAPACK_dsbev(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_ssbev LAPACK_GLOBAL(ssbev,SSBEV)
void LAPACK_ssbev(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_dsbev_2stage LAPACK_GLOBAL(dsbev_2stage,DSBEV_2STAGE)
void LAPACK_dsbev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssbev_2stage LAPACK_GLOBAL(ssbev_2stage,SSBEV_2STAGE)
void LAPACK_ssbev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsbevd LAPACK_GLOBAL(dsbevd,DSBEVD)
void LAPACK_dsbevd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssbevd LAPACK_GLOBAL(ssbevd,SSBEVD)
void LAPACK_ssbevd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsbevd_2stage LAPACK_GLOBAL(dsbevd_2stage,DSBEVD_2STAGE)
void LAPACK_dsbevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssbevd_2stage LAPACK_GLOBAL(ssbevd_2stage,SSBEVD_2STAGE)
void LAPACK_ssbevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsbevx LAPACK_GLOBAL(dsbevx,DSBEVX)
void LAPACK_dsbevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    double* Q, lapack_int const* ldq,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_ssbevx LAPACK_GLOBAL(ssbevx,SSBEVX)
void LAPACK_ssbevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    float* Q, lapack_int const* ldq,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_dsbevx_2stage LAPACK_GLOBAL(dsbevx_2stage,DSBEVX_2STAGE)
void LAPACK_dsbevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    double* Q, lapack_int const* ldq,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_ssbevx_2stage LAPACK_GLOBAL(ssbevx_2stage,SSBEVX_2STAGE)
void LAPACK_ssbevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    float* Q, lapack_int const* ldq,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_dsbgst LAPACK_GLOBAL(dsbgst,DSBGST)
void LAPACK_dsbgst(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    double* AB, lapack_int const* ldab,
    double const* BB, lapack_int const* ldbb,
    double* X, lapack_int const* ldx,
    double* work,
    lapack_int* info );

#define LAPACK_ssbgst LAPACK_GLOBAL(ssbgst,SSBGST)
void LAPACK_ssbgst(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    float* AB, lapack_int const* ldab,
    float const* BB, lapack_int const* ldbb,
    float* X, lapack_int const* ldx,
    float* work,
    lapack_int* info );

#define LAPACK_dsbgv LAPACK_GLOBAL(dsbgv,DSBGV)
void LAPACK_dsbgv(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    double* AB, lapack_int const* ldab,
    double* BB, lapack_int const* ldbb,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_ssbgv LAPACK_GLOBAL(ssbgv,SSBGV)
void LAPACK_ssbgv(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    float* AB, lapack_int const* ldab,
    float* BB, lapack_int const* ldbb,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_dsbgvd LAPACK_GLOBAL(dsbgvd,DSBGVD)
void LAPACK_dsbgvd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    double* AB, lapack_int const* ldab,
    double* BB, lapack_int const* ldbb,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssbgvd LAPACK_GLOBAL(ssbgvd,SSBGVD)
void LAPACK_ssbgvd(
    char const* jobz, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    float* AB, lapack_int const* ldab,
    float* BB, lapack_int const* ldbb,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsbgvx LAPACK_GLOBAL(dsbgvx,DSBGVX)
void LAPACK_dsbgvx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    double* AB, lapack_int const* ldab,
    double* BB, lapack_int const* ldbb,
    double* Q, lapack_int const* ldq,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_ssbgvx LAPACK_GLOBAL(ssbgvx,SSBGVX)
void LAPACK_ssbgvx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n, lapack_int const* ka, lapack_int const* kb,
    float* AB, lapack_int const* ldab,
    float* BB, lapack_int const* ldbb,
    float* Q, lapack_int const* ldq,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_dsbtrd LAPACK_GLOBAL(dsbtrd,DSBTRD)
void LAPACK_dsbtrd(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    double* AB, lapack_int const* ldab,
    double* D,
    double* E,
    double* Q, lapack_int const* ldq,
    double* work,
    lapack_int* info );

#define LAPACK_ssbtrd LAPACK_GLOBAL(ssbtrd,SSBTRD)
void LAPACK_ssbtrd(
    char const* vect, char const* uplo,
    lapack_int const* n, lapack_int const* kd,
    float* AB, lapack_int const* ldab,
    float* D,
    float* E,
    float* Q, lapack_int const* ldq,
    float* work,
    lapack_int* info );

#define LAPACK_dsfrk LAPACK_GLOBAL(dsfrk,DSFRK)
void LAPACK_dsfrk(
    char const* transr, char const* uplo, char const* trans,
    lapack_int const* n, lapack_int const* k,
    double const* alpha,
    double const* A, lapack_int const* lda,
    double const* beta,
    double* C );

#define LAPACK_ssfrk LAPACK_GLOBAL(ssfrk,SSFRK)
void LAPACK_ssfrk(
    char const* transr, char const* uplo, char const* trans,
    lapack_int const* n, lapack_int const* k,
    float const* alpha,
    float const* A, lapack_int const* lda,
    float const* beta,
    float* C );

#define LAPACK_cspcon LAPACK_GLOBAL(cspcon,CSPCON)
void LAPACK_cspcon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dspcon LAPACK_GLOBAL(dspcon,DSPCON)
void LAPACK_dspcon(
    char const* uplo,
    lapack_int const* n,
    double const* AP, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sspcon LAPACK_GLOBAL(sspcon,SSPCON)
void LAPACK_sspcon(
    char const* uplo,
    lapack_int const* n,
    float const* AP, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zspcon LAPACK_GLOBAL(zspcon,ZSPCON)
void LAPACK_zspcon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_dspev LAPACK_GLOBAL(dspev,DSPEV)
void LAPACK_dspev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    double* AP,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_sspev LAPACK_GLOBAL(sspev,SSPEV)
void LAPACK_sspev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    float* AP,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_dspevd LAPACK_GLOBAL(dspevd,DSPEVD)
void LAPACK_dspevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    double* AP,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_sspevd LAPACK_GLOBAL(sspevd,SSPEVD)
void LAPACK_sspevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    float* AP,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dspevx LAPACK_GLOBAL(dspevx,DSPEVX)
void LAPACK_dspevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    double* AP,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_sspevx LAPACK_GLOBAL(sspevx,SSPEVX)
void LAPACK_sspevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    float* AP,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_dspgst LAPACK_GLOBAL(dspgst,DSPGST)
void LAPACK_dspgst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    double* AP,
    double const* BP,
    lapack_int* info );

#define LAPACK_sspgst LAPACK_GLOBAL(sspgst,SSPGST)
void LAPACK_sspgst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    float* AP,
    float const* BP,
    lapack_int* info );

#define LAPACK_dspgv LAPACK_GLOBAL(dspgv,DSPGV)
void LAPACK_dspgv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    double* AP,
    double* BP,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_sspgv LAPACK_GLOBAL(sspgv,SSPGV)
void LAPACK_sspgv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    float* AP,
    float* BP,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_dspgvd LAPACK_GLOBAL(dspgvd,DSPGVD)
void LAPACK_dspgvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    double* AP,
    double* BP,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_sspgvd LAPACK_GLOBAL(sspgvd,SSPGVD)
void LAPACK_sspgvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    float* AP,
    float* BP,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dspgvx LAPACK_GLOBAL(dspgvx,DSPGVX)
void LAPACK_dspgvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    double* AP,
    double* BP,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_sspgvx LAPACK_GLOBAL(sspgvx,SSPGVX)
void LAPACK_sspgvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    float* AP,
    float* BP,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_csprfs LAPACK_GLOBAL(csprfs,CSPRFS)
void LAPACK_csprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float const* AFP, lapack_int const* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dsprfs LAPACK_GLOBAL(dsprfs,DSPRFS)
void LAPACK_dsprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* AP,
    double const* AFP, lapack_int const* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ssprfs LAPACK_GLOBAL(ssprfs,SSPRFS)
void LAPACK_ssprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* AP,
    float const* AFP, lapack_int const* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zsprfs LAPACK_GLOBAL(zsprfs,ZSPRFS)
void LAPACK_zsprfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double const* AFP, lapack_int const* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_cspsv LAPACK_GLOBAL(cspsv,CSPSV)
void LAPACK_cspsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* AP, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dspsv LAPACK_GLOBAL(dspsv,DSPSV)
void LAPACK_dspsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* AP, lapack_int* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_sspsv LAPACK_GLOBAL(sspsv,SSPSV)
void LAPACK_sspsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* AP, lapack_int* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zspsv LAPACK_GLOBAL(zspsv,ZSPSV)
void LAPACK_zspsv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* AP, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_cspsvx LAPACK_GLOBAL(cspsvx,CSPSVX)
void LAPACK_cspsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float* AFP, lapack_int* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dspsvx LAPACK_GLOBAL(dspsvx,DSPSVX)
void LAPACK_dspsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* AP,
    double* AFP, lapack_int* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sspsvx LAPACK_GLOBAL(sspsvx,SSPSVX)
void LAPACK_sspsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* AP,
    float* AFP, lapack_int* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zspsvx LAPACK_GLOBAL(zspsvx,ZSPSVX)
void LAPACK_zspsvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double* AFP, lapack_int* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_dsptrd LAPACK_GLOBAL(dsptrd,DSPTRD)
void LAPACK_dsptrd(
    char const* uplo,
    lapack_int const* n,
    double* AP,
    double* D,
    double* E,
    double* tau,
    lapack_int* info );

#define LAPACK_ssptrd LAPACK_GLOBAL(ssptrd,SSPTRD)
void LAPACK_ssptrd(
    char const* uplo,
    lapack_int const* n,
    float* AP,
    float* D,
    float* E,
    float* tau,
    lapack_int* info );

#define LAPACK_csptrf LAPACK_GLOBAL(csptrf,CSPTRF)
void LAPACK_csptrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_dsptrf LAPACK_GLOBAL(dsptrf,DSPTRF)
void LAPACK_dsptrf(
    char const* uplo,
    lapack_int const* n,
    double* AP, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_ssptrf LAPACK_GLOBAL(ssptrf,SSPTRF)
void LAPACK_ssptrf(
    char const* uplo,
    lapack_int const* n,
    float* AP, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_zsptrf LAPACK_GLOBAL(zsptrf,ZSPTRF)
void LAPACK_zsptrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP, lapack_int* ipiv,
    lapack_int* info );

#define LAPACK_csptri LAPACK_GLOBAL(csptri,CSPTRI)
void LAPACK_csptri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* AP, lapack_int const* ipiv,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dsptri LAPACK_GLOBAL(dsptri,DSPTRI)
void LAPACK_dsptri(
    char const* uplo,
    lapack_int const* n,
    double* AP, lapack_int const* ipiv,
    double* work,
    lapack_int* info );

#define LAPACK_ssptri LAPACK_GLOBAL(ssptri,SSPTRI)
void LAPACK_ssptri(
    char const* uplo,
    lapack_int const* n,
    float* AP, lapack_int const* ipiv,
    float* work,
    lapack_int* info );

#define LAPACK_zsptri LAPACK_GLOBAL(zsptri,ZSPTRI)
void LAPACK_zsptri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP, lapack_int const* ipiv,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_csptrs LAPACK_GLOBAL(csptrs,CSPTRS)
void LAPACK_csptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsptrs LAPACK_GLOBAL(dsptrs,DSPTRS)
void LAPACK_dsptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* AP, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ssptrs LAPACK_GLOBAL(ssptrs,SSPTRS)
void LAPACK_ssptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* AP, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zsptrs LAPACK_GLOBAL(zsptrs,ZSPTRS)
void LAPACK_zsptrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dstebz LAPACK_GLOBAL(dstebz,DSTEBZ)
void LAPACK_dstebz(
    char const* range, char const* order,
    lapack_int const* n,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol,
    double const* D,
    double const* E, lapack_int* m, lapack_int* nsplit,
    double* W, lapack_int* IBLOCK, lapack_int* ISPLIT,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_sstebz LAPACK_GLOBAL(sstebz,SSTEBZ)
void LAPACK_sstebz(
    char const* range, char const* order,
    lapack_int const* n,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol,
    float const* D,
    float const* E, lapack_int* m, lapack_int* nsplit,
    float* W, lapack_int* IBLOCK, lapack_int* ISPLIT,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cstedc LAPACK_GLOBAL(cstedc,CSTEDC)
void LAPACK_cstedc(
    char const* compz,
    lapack_int const* n,
    float* D,
    float* E,
    lapack_complex_float* Z, lapack_int const* ldz,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dstedc LAPACK_GLOBAL(dstedc,DSTEDC)
void LAPACK_dstedc(
    char const* compz,
    lapack_int const* n,
    double* D,
    double* E,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_sstedc LAPACK_GLOBAL(sstedc,SSTEDC)
void LAPACK_sstedc(
    char const* compz,
    lapack_int const* n,
    float* D,
    float* E,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zstedc LAPACK_GLOBAL(zstedc,ZSTEDC)
void LAPACK_zstedc(
    char const* compz,
    lapack_int const* n,
    double* D,
    double* E,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_cstegr LAPACK_GLOBAL(cstegr,CSTEGR)
void LAPACK_cstegr(
    char const* jobz, char const* range,
    lapack_int const* n,
    float* D,
    float* E,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dstegr LAPACK_GLOBAL(dstegr,DSTEGR)
void LAPACK_dstegr(
    char const* jobz, char const* range,
    lapack_int const* n,
    double* D,
    double* E,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_sstegr LAPACK_GLOBAL(sstegr,SSTEGR)
void LAPACK_sstegr(
    char const* jobz, char const* range,
    lapack_int const* n,
    float* D,
    float* E,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zstegr LAPACK_GLOBAL(zstegr,ZSTEGR)
void LAPACK_zstegr(
    char const* jobz, char const* range,
    lapack_int const* n,
    double* D,
    double* E,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_cstein LAPACK_GLOBAL(cstein,CSTEIN)
void LAPACK_cstein(
    lapack_int const* n,
    float const* D,
    float const* E, lapack_int const* m,
    float const* W, lapack_int const* IBLOCK, lapack_int const* ISPLIT,
    lapack_complex_float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_dstein LAPACK_GLOBAL(dstein,DSTEIN)
void LAPACK_dstein(
    lapack_int const* n,
    double const* D,
    double const* E, lapack_int const* m,
    double const* W, lapack_int const* IBLOCK, lapack_int const* ISPLIT,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_sstein LAPACK_GLOBAL(sstein,SSTEIN)
void LAPACK_sstein(
    lapack_int const* n,
    float const* D,
    float const* E, lapack_int const* m,
    float const* W, lapack_int const* IBLOCK, lapack_int const* ISPLIT,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_zstein LAPACK_GLOBAL(zstein,ZSTEIN)
void LAPACK_zstein(
    lapack_int const* n,
    double const* D,
    double const* E, lapack_int const* m,
    double const* W, lapack_int const* IBLOCK, lapack_int const* ISPLIT,
    lapack_complex_double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_cstemr LAPACK_GLOBAL(cstemr,CSTEMR)
void LAPACK_cstemr(
    char const* jobz, char const* range,
    lapack_int const* n,
    float* D,
    float* E,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* m,
    float* W,
    lapack_complex_float* Z, lapack_int const* ldz, lapack_int const* nzc, lapack_int* ISUPPZ, lapack_logical* tryrac,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dstemr LAPACK_GLOBAL(dstemr,DSTEMR)
void LAPACK_dstemr(
    char const* jobz, char const* range,
    lapack_int const* n,
    double* D,
    double* E,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz, lapack_int const* nzc, lapack_int* ISUPPZ, lapack_logical* tryrac,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_sstemr LAPACK_GLOBAL(sstemr,SSTEMR)
void LAPACK_sstemr(
    char const* jobz, char const* range,
    lapack_int const* n,
    float* D,
    float* E,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz, lapack_int const* nzc, lapack_int* ISUPPZ, lapack_logical* tryrac,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_zstemr LAPACK_GLOBAL(zstemr,ZSTEMR)
void LAPACK_zstemr(
    char const* jobz, char const* range,
    lapack_int const* n,
    double* D,
    double* E,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz, lapack_int const* nzc, lapack_int* ISUPPZ, lapack_logical* tryrac,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_csteqr LAPACK_GLOBAL(csteqr,CSTEQR)
void LAPACK_csteqr(
    char const* compz,
    lapack_int const* n,
    float* D,
    float* E,
    lapack_complex_float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_dsteqr LAPACK_GLOBAL(dsteqr,DSTEQR)
void LAPACK_dsteqr(
    char const* compz,
    lapack_int const* n,
    double* D,
    double* E,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_ssteqr LAPACK_GLOBAL(ssteqr,SSTEQR)
void LAPACK_ssteqr(
    char const* compz,
    lapack_int const* n,
    float* D,
    float* E,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_zsteqr LAPACK_GLOBAL(zsteqr,ZSTEQR)
void LAPACK_zsteqr(
    char const* compz,
    lapack_int const* n,
    double* D,
    double* E,
    lapack_complex_double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_dsterf LAPACK_GLOBAL(dsterf,DSTERF)
void LAPACK_dsterf(
    lapack_int const* n,
    double* D,
    double* E,
    lapack_int* info );

#define LAPACK_ssterf LAPACK_GLOBAL(ssterf,SSTERF)
void LAPACK_ssterf(
    lapack_int const* n,
    float* D,
    float* E,
    lapack_int* info );

#define LAPACK_dstev LAPACK_GLOBAL(dstev,DSTEV)
void LAPACK_dstev(
    char const* jobz,
    lapack_int const* n,
    double* D,
    double* E,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* info );

#define LAPACK_sstev LAPACK_GLOBAL(sstev,SSTEV)
void LAPACK_sstev(
    char const* jobz,
    lapack_int const* n,
    float* D,
    float* E,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* info );

#define LAPACK_dstevd LAPACK_GLOBAL(dstevd,DSTEVD)
void LAPACK_dstevd(
    char const* jobz,
    lapack_int const* n,
    double* D,
    double* E,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_sstevd LAPACK_GLOBAL(sstevd,SSTEVD)
void LAPACK_sstevd(
    char const* jobz,
    lapack_int const* n,
    float* D,
    float* E,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dstevr LAPACK_GLOBAL(dstevr,DSTEVR)
void LAPACK_dstevr(
    char const* jobz, char const* range,
    lapack_int const* n,
    double* D,
    double* E,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_sstevr LAPACK_GLOBAL(sstevr,SSTEVR)
void LAPACK_sstevr(
    char const* jobz, char const* range,
    lapack_int const* n,
    float* D,
    float* E,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dstevx LAPACK_GLOBAL(dstevx,DSTEVX)
void LAPACK_dstevx(
    char const* jobz, char const* range,
    lapack_int const* n,
    double* D,
    double* E,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_sstevx LAPACK_GLOBAL(sstevx,SSTEVX)
void LAPACK_sstevx(
    char const* jobz, char const* range,
    lapack_int const* n,
    float* D,
    float* E,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_csycon LAPACK_GLOBAL(csycon,CSYCON)
void LAPACK_csycon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dsycon LAPACK_GLOBAL(dsycon,DSYCON)
void LAPACK_dsycon(
    char const* uplo,
    lapack_int const* n,
    double const* A, lapack_int const* lda, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ssycon LAPACK_GLOBAL(ssycon,SSYCON)
void LAPACK_ssycon(
    char const* uplo,
    lapack_int const* n,
    float const* A, lapack_int const* lda, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zsycon LAPACK_GLOBAL(zsycon,ZSYCON)
void LAPACK_zsycon(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_csycon_3 LAPACK_GLOBAL(csycon_3,CSYCON_3)
void LAPACK_csycon_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* E, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dsycon_3 LAPACK_GLOBAL(dsycon_3,DSYCON_3)
void LAPACK_dsycon_3(
    char const* uplo,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double const* E, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ssycon_3 LAPACK_GLOBAL(ssycon_3,SSYCON_3)
void LAPACK_ssycon_3(
    char const* uplo,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float const* E, lapack_int const* ipiv,
    float const* anorm,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zsycon_3 LAPACK_GLOBAL(zsycon_3,ZSYCON_3)
void LAPACK_zsycon_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* E, lapack_int const* ipiv,
    double const* anorm,
    double* rcond,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_csyconv LAPACK_GLOBAL(csyconv,CSYCONV)
void LAPACK_csyconv(
    char const* uplo, char const* way,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* E,
    lapack_int* info );

#define LAPACK_dsyconv LAPACK_GLOBAL(dsyconv,DSYCONV)
void LAPACK_dsyconv(
    char const* uplo, char const* way,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int const* ipiv,
    double* E,
    lapack_int* info );

#define LAPACK_ssyconv LAPACK_GLOBAL(ssyconv,SSYCONV)
void LAPACK_ssyconv(
    char const* uplo, char const* way,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int const* ipiv,
    float* E,
    lapack_int* info );

#define LAPACK_zsyconv LAPACK_GLOBAL(zsyconv,ZSYCONV)
void LAPACK_zsyconv(
    char const* uplo, char const* way,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* E,
    lapack_int* info );

#define LAPACK_csyequb LAPACK_GLOBAL(csyequb,CSYEQUB)
void LAPACK_csyequb(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* S,
    float* scond,
    float* amax,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dsyequb LAPACK_GLOBAL(dsyequb,DSYEQUB)
void LAPACK_dsyequb(
    char const* uplo,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* S,
    double* scond,
    double* amax,
    double* work,
    lapack_int* info );

#define LAPACK_ssyequb LAPACK_GLOBAL(ssyequb,SSYEQUB)
void LAPACK_ssyequb(
    char const* uplo,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* S,
    float* scond,
    float* amax,
    float* work,
    lapack_int* info );

#define LAPACK_zsyequb LAPACK_GLOBAL(zsyequb,ZSYEQUB)
void LAPACK_zsyequb(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* S,
    double* scond,
    double* amax,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_dsyev LAPACK_GLOBAL(dsyev,DSYEV)
void LAPACK_dsyev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* W,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssyev LAPACK_GLOBAL(ssyev,SSYEV)
void LAPACK_ssyev(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* W,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsyev_2stage LAPACK_GLOBAL(dsyev_2stage,DSYEV_2STAGE)
void LAPACK_dsyev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* W,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssyev_2stage LAPACK_GLOBAL(ssyev_2stage,SSYEV_2STAGE)
void LAPACK_ssyev_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* W,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsyevd LAPACK_GLOBAL(dsyevd,DSYEVD)
void LAPACK_dsyevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* W,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssyevd LAPACK_GLOBAL(ssyevd,SSYEVD)
void LAPACK_ssyevd(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* W,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsyevd_2stage LAPACK_GLOBAL(dsyevd_2stage,DSYEVD_2STAGE)
void LAPACK_dsyevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* W,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssyevd_2stage LAPACK_GLOBAL(ssyevd_2stage,SSYEVD_2STAGE)
void LAPACK_ssyevd_2stage(
    char const* jobz, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* W,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsyevr LAPACK_GLOBAL(dsyevr,DSYEVR)
void LAPACK_dsyevr(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssyevr LAPACK_GLOBAL(ssyevr,SSYEVR)
void LAPACK_ssyevr(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsyevr_2stage LAPACK_GLOBAL(dsyevr_2stage,DSYEVR_2STAGE)
void LAPACK_dsyevr_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssyevr_2stage LAPACK_GLOBAL(ssyevr_2stage,SSYEVR_2STAGE)
void LAPACK_ssyevr_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz, lapack_int* ISUPPZ,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsyevx LAPACK_GLOBAL(dsyevx,DSYEVX)
void LAPACK_dsyevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_ssyevx LAPACK_GLOBAL(ssyevx,SSYEVX)
void LAPACK_ssyevx(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_dsyevx_2stage LAPACK_GLOBAL(dsyevx_2stage,DSYEVX_2STAGE)
void LAPACK_dsyevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_ssyevx_2stage LAPACK_GLOBAL(ssyevx_2stage,SSYEVX_2STAGE)
void LAPACK_ssyevx_2stage(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_dsygst LAPACK_GLOBAL(dsygst,DSYGST)
void LAPACK_dsygst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ssygst LAPACK_GLOBAL(ssygst,SSYGST)
void LAPACK_ssygst(
    lapack_int const* itype, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsygv LAPACK_GLOBAL(dsygv,DSYGV)
void LAPACK_dsygv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* W,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssygv LAPACK_GLOBAL(ssygv,SSYGV)
void LAPACK_ssygv(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* W,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsygv_2stage LAPACK_GLOBAL(dsygv_2stage,DSYGV_2STAGE)
void LAPACK_dsygv_2stage(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* W,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssygv_2stage LAPACK_GLOBAL(ssygv_2stage,SSYGV_2STAGE)
void LAPACK_ssygv_2stage(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* W,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsygvd LAPACK_GLOBAL(dsygvd,DSYGVD)
void LAPACK_dsygvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* W,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ssygvd LAPACK_GLOBAL(ssygvd,SSYGVD)
void LAPACK_ssygvd(
    lapack_int const* itype, char const* jobz, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* W,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dsygvx LAPACK_GLOBAL(dsygvx,DSYGVX)
void LAPACK_dsygvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    double* Z, lapack_int const* ldz,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_ssygvx LAPACK_GLOBAL(ssygvx,SSYGVX)
void LAPACK_ssygvx(
    lapack_int const* itype, char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float const* vl,
    float const* vu, lapack_int const* il, lapack_int const* iu,
    float const* abstol, lapack_int* m,
    float* W,
    float* Z, lapack_int const* ldz,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info );

#define LAPACK_csyr LAPACK_GLOBAL(csyr,CSYR)
void LAPACK_csyr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* alpha,
    lapack_complex_float const* X, lapack_int const* incx,
    lapack_complex_float* A, lapack_int const* lda );

#define LAPACK_zsyr LAPACK_GLOBAL(zsyr,ZSYR)
void LAPACK_zsyr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* alpha,
    lapack_complex_double const* X, lapack_int const* incx,
    lapack_complex_double* A, lapack_int const* lda );

#define LAPACK_csyrfs LAPACK_GLOBAL(csyrfs,CSYRFS)
void LAPACK_csyrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dsyrfs LAPACK_GLOBAL(dsyrfs,DSYRFS)
void LAPACK_dsyrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ssyrfs LAPACK_GLOBAL(ssyrfs,SSYRFS)
void LAPACK_ssyrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zsyrfs LAPACK_GLOBAL(zsyrfs,ZSYRFS)
void LAPACK_zsyrfs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_csyrfsx LAPACK_GLOBAL(csyrfsx,CSYRFSX)
void LAPACK_csyrfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    float* S,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dsyrfsx LAPACK_GLOBAL(dsyrfsx,DSYRFSX)
void LAPACK_dsyrfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    double* S,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ssyrfsx LAPACK_GLOBAL(ssyrfsx,SSYRFSX)
void LAPACK_ssyrfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    float* S,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zsyrfsx LAPACK_GLOBAL(zsyrfsx,ZSYRFSX)
void LAPACK_zsyrfsx(
    char const* uplo, char const* equed,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* AF, lapack_int const* ldaf, lapack_int const* ipiv,
    double* S,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_csysv LAPACK_GLOBAL(csysv,CSYSV)
void LAPACK_csysv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsysv LAPACK_GLOBAL(dsysv,DSYSV)
void LAPACK_dsysv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssysv LAPACK_GLOBAL(ssysv,SSYSV)
void LAPACK_ssysv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsysv LAPACK_GLOBAL(zsysv,ZSYSV)
void LAPACK_zsysv(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csysv_aa LAPACK_GLOBAL(csysv_aa,CSYSV_AA)
void LAPACK_csysv_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsysv_aa LAPACK_GLOBAL(dsysv_aa,DSYSV_AA)
void LAPACK_dsysv_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssysv_aa LAPACK_GLOBAL(ssysv_aa,SSYSV_AA)
void LAPACK_ssysv_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsysv_aa LAPACK_GLOBAL(zsysv_aa,ZSYSV_AA)
void LAPACK_zsysv_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csysv_aa_2stage LAPACK_GLOBAL(csysv_aa_2stage,CSYSV_AA_2STAGE)
void LAPACK_csysv_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsysv_aa_2stage LAPACK_GLOBAL(dsysv_aa_2stage,DSYSV_AA_2STAGE)
void LAPACK_dsysv_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssysv_aa_2stage LAPACK_GLOBAL(ssysv_aa_2stage,SSYSV_AA_2STAGE)
void LAPACK_ssysv_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsysv_aa_2stage LAPACK_GLOBAL(zsysv_aa_2stage,ZSYSV_AA_2STAGE)
void LAPACK_zsysv_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csysv_rk LAPACK_GLOBAL(csysv_rk,CSYSV_RK)
void LAPACK_csysv_rk(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* E, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsysv_rk LAPACK_GLOBAL(dsysv_rk,DSYSV_RK)
void LAPACK_dsysv_rk(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* E, lapack_int* ipiv,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssysv_rk LAPACK_GLOBAL(ssysv_rk,SSYSV_RK)
void LAPACK_ssysv_rk(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* E, lapack_int* ipiv,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsysv_rk LAPACK_GLOBAL(zsysv_rk,ZSYSV_RK)
void LAPACK_zsysv_rk(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* E, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csysv_rook LAPACK_GLOBAL(csysv_rook,CSYSV_ROOK)
void LAPACK_csysv_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsysv_rook LAPACK_GLOBAL(dsysv_rook,DSYSV_ROOK)
void LAPACK_dsysv_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssysv_rook LAPACK_GLOBAL(ssysv_rook,SSYSV_ROOK)
void LAPACK_ssysv_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsysv_rook LAPACK_GLOBAL(zsysv_rook,ZSYSV_ROOK)
void LAPACK_zsysv_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csysvx LAPACK_GLOBAL(csysvx,CSYSVX)
void LAPACK_csysvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, lapack_int* ipiv,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dsysvx LAPACK_GLOBAL(dsysvx,DSYSVX)
void LAPACK_dsysvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double* AF, lapack_int const* ldaf, lapack_int* ipiv,
    double const* B, lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ssysvx LAPACK_GLOBAL(ssysvx,SSYSVX)
void LAPACK_ssysvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float* AF, lapack_int const* ldaf, lapack_int* ipiv,
    float const* B, lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* ferr,
    float* berr,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zsysvx LAPACK_GLOBAL(zsysvx,ZSYSVX)
void LAPACK_zsysvx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, lapack_int* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_csysvxx LAPACK_GLOBAL(csysvxx,CSYSVXX)
void LAPACK_csysvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    float* S,
    lapack_complex_float* B,
    lapack_int const* ldb,
    lapack_complex_float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dsysvxx LAPACK_GLOBAL(dsysvxx,DSYSVXX)
void LAPACK_dsysvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda,
    double* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    double* S,
    double* B,
    lapack_int const* ldb,
    double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ssysvxx LAPACK_GLOBAL(ssysvxx,SSYSVXX)
void LAPACK_ssysvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda,
    float* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    float* S,
    float* B,
    lapack_int const* ldb,
    float* X, lapack_int const* ldx,
    float* rcond,
    float* rpvgrw,
    float* berr, lapack_int const* n_err_bnds,
    float* err_bnds_norm,
    float* err_bnds_comp, lapack_int const* nparams,
    float* params,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zsysvxx LAPACK_GLOBAL(zsysvxx,ZSYSVXX)
void LAPACK_zsysvxx(
    char const* fact, char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* AF, lapack_int const* ldaf, lapack_int* ipiv, char* equed,
    double* S,
    lapack_complex_double* B,
    lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* rpvgrw,
    double* berr, lapack_int const* n_err_bnds,
    double* err_bnds_norm,
    double* err_bnds_comp, lapack_int const* nparams,
    double* params,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_csyswapr LAPACK_GLOBAL(csyswapr,CSYSWAPR)
void LAPACK_csyswapr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* i1, lapack_int const* i2 );

#define LAPACK_dsyswapr LAPACK_GLOBAL(dsyswapr,DSYSWAPR)
void LAPACK_dsyswapr(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int const* i1, lapack_int const* i2 );

#define LAPACK_ssyswapr LAPACK_GLOBAL(ssyswapr,SSYSWAPR)
void LAPACK_ssyswapr(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int const* i1, lapack_int const* i2 );

#define LAPACK_zsyswapr LAPACK_GLOBAL(zsyswapr,ZSYSWAPR)
void LAPACK_zsyswapr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* i1, lapack_int const* i2 );

#define LAPACK_dsytrd LAPACK_GLOBAL(dsytrd,DSYTRD)
void LAPACK_dsytrd(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* D,
    double* E,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrd LAPACK_GLOBAL(ssytrd,SSYTRD)
void LAPACK_ssytrd(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* D,
    float* E,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytrd_2stage LAPACK_GLOBAL(dsytrd_2stage,DSYTRD_2STAGE)
void LAPACK_dsytrd_2stage(
    char const* vect, char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* D,
    double* E,
    double* tau,
    double* HOUS2, lapack_int const* lhous2,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrd_2stage LAPACK_GLOBAL(ssytrd_2stage,SSYTRD_2STAGE)
void LAPACK_ssytrd_2stage(
    char const* vect, char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* D,
    float* E,
    float* tau,
    float* HOUS2, lapack_int const* lhous2,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytrf LAPACK_GLOBAL(csytrf,CSYTRF)
void LAPACK_csytrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytrf LAPACK_GLOBAL(dsytrf,DSYTRF)
void LAPACK_dsytrf(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrf LAPACK_GLOBAL(ssytrf,SSYTRF)
void LAPACK_ssytrf(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytrf LAPACK_GLOBAL(zsytrf,ZSYTRF)
void LAPACK_zsytrf(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytrf_aa LAPACK_GLOBAL(csytrf_aa,CSYTRF_AA)
void LAPACK_csytrf_aa(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytrf_aa LAPACK_GLOBAL(dsytrf_aa,DSYTRF_AA)
void LAPACK_dsytrf_aa(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrf_aa LAPACK_GLOBAL(ssytrf_aa,SSYTRF_AA)
void LAPACK_ssytrf_aa(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytrf_aa LAPACK_GLOBAL(zsytrf_aa,ZSYTRF_AA)
void LAPACK_zsytrf_aa(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytrf_aa_2stage LAPACK_GLOBAL(csytrf_aa_2stage,CSYTRF_AA_2STAGE)
void LAPACK_csytrf_aa_2stage(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytrf_aa_2stage LAPACK_GLOBAL(dsytrf_aa_2stage,DSYTRF_AA_2STAGE)
void LAPACK_dsytrf_aa_2stage(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrf_aa_2stage LAPACK_GLOBAL(ssytrf_aa_2stage,SSYTRF_AA_2STAGE)
void LAPACK_ssytrf_aa_2stage(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytrf_aa_2stage LAPACK_GLOBAL(zsytrf_aa_2stage,ZSYTRF_AA_2STAGE)
void LAPACK_zsytrf_aa_2stage(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* TB, lapack_int const* ltb, lapack_int* ipiv, lapack_int* ipiv2,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytrf_rk LAPACK_GLOBAL(csytrf_rk,CSYTRF_RK)
void LAPACK_csytrf_rk(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* E, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytrf_rk LAPACK_GLOBAL(dsytrf_rk,DSYTRF_RK)
void LAPACK_dsytrf_rk(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double* E, lapack_int* ipiv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrf_rk LAPACK_GLOBAL(ssytrf_rk,SSYTRF_RK)
void LAPACK_ssytrf_rk(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float* E, lapack_int* ipiv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytrf_rk LAPACK_GLOBAL(zsytrf_rk,ZSYTRF_RK)
void LAPACK_zsytrf_rk(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* E, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytrf_rook LAPACK_GLOBAL(csytrf_rook,CSYTRF_ROOK)
void LAPACK_csytrf_rook(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytrf_rook LAPACK_GLOBAL(dsytrf_rook,DSYTRF_ROOK)
void LAPACK_dsytrf_rook(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int* ipiv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrf_rook LAPACK_GLOBAL(ssytrf_rook,SSYTRF_ROOK)
void LAPACK_ssytrf_rook(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int* ipiv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytrf_rook LAPACK_GLOBAL(zsytrf_rook,ZSYTRF_ROOK)
void LAPACK_zsytrf_rook(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytri LAPACK_GLOBAL(csytri,CSYTRI)
void LAPACK_csytri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dsytri LAPACK_GLOBAL(dsytri,DSYTRI)
void LAPACK_dsytri(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int const* ipiv,
    double* work,
    lapack_int* info );

#define LAPACK_ssytri LAPACK_GLOBAL(ssytri,SSYTRI)
void LAPACK_ssytri(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int const* ipiv,
    float* work,
    lapack_int* info );

#define LAPACK_zsytri LAPACK_GLOBAL(zsytri,ZSYTRI)
void LAPACK_zsytri(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_csytri2 LAPACK_GLOBAL(csytri2,CSYTRI2)
void LAPACK_csytri2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytri2 LAPACK_GLOBAL(dsytri2,DSYTRI2)
void LAPACK_dsytri2(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int const* ipiv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytri2 LAPACK_GLOBAL(ssytri2,SSYTRI2)
void LAPACK_ssytri2(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int const* ipiv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytri2 LAPACK_GLOBAL(zsytri2,ZSYTRI2)
void LAPACK_zsytri2(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytri2x LAPACK_GLOBAL(csytri2x,CSYTRI2X)
void LAPACK_csytri2x(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* work, lapack_int const* nb,
    lapack_int* info );

#define LAPACK_dsytri2x LAPACK_GLOBAL(dsytri2x,DSYTRI2X)
void LAPACK_dsytri2x(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda, lapack_int const* ipiv,
    double* work, lapack_int const* nb,
    lapack_int* info );

#define LAPACK_ssytri2x LAPACK_GLOBAL(ssytri2x,SSYTRI2X)
void LAPACK_ssytri2x(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda, lapack_int const* ipiv,
    float* work, lapack_int const* nb,
    lapack_int* info );

#define LAPACK_zsytri2x LAPACK_GLOBAL(zsytri2x,ZSYTRI2X)
void LAPACK_zsytri2x(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* work, lapack_int const* nb,
    lapack_int* info );

#define LAPACK_csytri_3 LAPACK_GLOBAL(csytri_3,CSYTRI_3)
void LAPACK_csytri_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* E, lapack_int const* ipiv,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytri_3 LAPACK_GLOBAL(dsytri_3,DSYTRI_3)
void LAPACK_dsytri_3(
    char const* uplo,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    double const* E, lapack_int const* ipiv,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytri_3 LAPACK_GLOBAL(ssytri_3,SSYTRI_3)
void LAPACK_ssytri_3(
    char const* uplo,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    float const* E, lapack_int const* ipiv,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytri_3 LAPACK_GLOBAL(zsytri_3,ZSYTRI_3)
void LAPACK_zsytri_3(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* E, lapack_int const* ipiv,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytrs LAPACK_GLOBAL(csytrs,CSYTRS)
void LAPACK_csytrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsytrs LAPACK_GLOBAL(dsytrs,DSYTRS)
void LAPACK_dsytrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ssytrs LAPACK_GLOBAL(ssytrs,SSYTRS)
void LAPACK_ssytrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zsytrs LAPACK_GLOBAL(zsytrs,ZSYTRS)
void LAPACK_zsytrs(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_csytrs2 LAPACK_GLOBAL(csytrs2,CSYTRS2)
void LAPACK_csytrs2(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dsytrs2 LAPACK_GLOBAL(dsytrs2,DSYTRS2)
void LAPACK_dsytrs2(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double* A, lapack_int const* lda, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    double* work,
    lapack_int* info );

#define LAPACK_ssytrs2 LAPACK_GLOBAL(ssytrs2,SSYTRS2)
void LAPACK_ssytrs2(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float* A, lapack_int const* lda, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    float* work,
    lapack_int* info );

#define LAPACK_zsytrs2 LAPACK_GLOBAL(zsytrs2,ZSYTRS2)
void LAPACK_zsytrs2(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_csytrs_3 LAPACK_GLOBAL(csytrs_3,CSYTRS_3)
void LAPACK_csytrs_3(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* E, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsytrs_3 LAPACK_GLOBAL(dsytrs_3,DSYTRS_3)
void LAPACK_dsytrs_3(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* E, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ssytrs_3 LAPACK_GLOBAL(ssytrs_3,SSYTRS_3)
void LAPACK_ssytrs_3(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* E, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zsytrs_3 LAPACK_GLOBAL(zsytrs_3,ZSYTRS_3)
void LAPACK_zsytrs_3(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* E, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_csytrs_aa LAPACK_GLOBAL(csytrs_aa,CSYTRS_AA)
void LAPACK_csytrs_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dsytrs_aa LAPACK_GLOBAL(dsytrs_aa,DSYTRS_AA)
void LAPACK_dsytrs_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ssytrs_aa LAPACK_GLOBAL(ssytrs_aa,SSYTRS_AA)
void LAPACK_ssytrs_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zsytrs_aa LAPACK_GLOBAL(zsytrs_aa,ZSYTRS_AA)
void LAPACK_zsytrs_aa(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_csytrs_aa_2stage LAPACK_GLOBAL(csytrs_aa_2stage,CSYTRS_AA_2STAGE)
void LAPACK_csytrs_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* TB, lapack_int const* ltb, lapack_int const* ipiv, lapack_int const* ipiv2,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsytrs_aa_2stage LAPACK_GLOBAL(dsytrs_aa_2stage,DSYTRS_AA_2STAGE)
void LAPACK_dsytrs_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double* TB, lapack_int const* ltb, lapack_int const* ipiv, lapack_int const* ipiv2,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ssytrs_aa_2stage LAPACK_GLOBAL(ssytrs_aa_2stage,SSYTRS_AA_2STAGE)
void LAPACK_ssytrs_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float* TB, lapack_int const* ltb, lapack_int const* ipiv, lapack_int const* ipiv2,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zsytrs_aa_2stage LAPACK_GLOBAL(zsytrs_aa_2stage,ZSYTRS_AA_2STAGE)
void LAPACK_zsytrs_aa_2stage(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* TB, lapack_int const* ltb, lapack_int const* ipiv, lapack_int const* ipiv2,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_csytrs_rook LAPACK_GLOBAL(csytrs_rook,CSYTRS_ROOK)
void LAPACK_csytrs_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dsytrs_rook LAPACK_GLOBAL(dsytrs_rook,DSYTRS_ROOK)
void LAPACK_dsytrs_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ssytrs_rook LAPACK_GLOBAL(ssytrs_rook,SSYTRS_ROOK)
void LAPACK_ssytrs_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda, lapack_int const* ipiv,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zsytrs_rook LAPACK_GLOBAL(zsytrs_rook,ZSYTRS_ROOK)
void LAPACK_zsytrs_rook(
    char const* uplo,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda, lapack_int const* ipiv,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ctbcon LAPACK_GLOBAL(ctbcon,CTBCON)
void LAPACK_ctbcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_float const* AB, lapack_int const* ldab,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtbcon LAPACK_GLOBAL(dtbcon,DTBCON)
void LAPACK_dtbcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* kd,
    double const* AB, lapack_int const* ldab,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_stbcon LAPACK_GLOBAL(stbcon,STBCON)
void LAPACK_stbcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* kd,
    float const* AB, lapack_int const* ldab,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztbcon LAPACK_GLOBAL(ztbcon,ZTBCON)
void LAPACK_ztbcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n, lapack_int const* kd,
    lapack_complex_double const* AB, lapack_int const* ldab,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctbrfs LAPACK_GLOBAL(ctbrfs,CTBRFS)
void LAPACK_ctbrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_float const* AB, lapack_int const* ldab,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float const* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtbrfs LAPACK_GLOBAL(dtbrfs,DTBRFS)
void LAPACK_dtbrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab,
    double const* B, lapack_int const* ldb,
    double const* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_stbrfs LAPACK_GLOBAL(stbrfs,STBRFS)
void LAPACK_stbrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    float const* AB, lapack_int const* ldab,
    float const* B, lapack_int const* ldb,
    float const* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztbrfs LAPACK_GLOBAL(ztbrfs,ZTBRFS)
void LAPACK_ztbrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_double const* AB, lapack_int const* ldab,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double const* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctbtrs LAPACK_GLOBAL(ctbtrs,CTBTRS)
void LAPACK_ctbtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_float const* AB, lapack_int const* ldab,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dtbtrs LAPACK_GLOBAL(dtbtrs,DTBTRS)
void LAPACK_dtbtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_stbtrs LAPACK_GLOBAL(stbtrs,STBTRS)
void LAPACK_stbtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    float const* AB, lapack_int const* ldab,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ztbtrs LAPACK_GLOBAL(ztbtrs,ZTBTRS)
void LAPACK_ztbtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* kd, lapack_int const* nrhs,
    lapack_complex_double const* AB, lapack_int const* ldab,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ctfsm LAPACK_GLOBAL(ctfsm,CTFSM)
void LAPACK_ctfsm(
    char const* transr, char const* side, char const* uplo, char const* trans, char const* diag,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* alpha,
    lapack_complex_float const* A,
    lapack_complex_float* B, lapack_int const* ldb );

#define LAPACK_dtfsm LAPACK_GLOBAL(dtfsm,DTFSM)
void LAPACK_dtfsm(
    char const* transr, char const* side, char const* uplo, char const* trans, char const* diag,
    lapack_int const* m, lapack_int const* n,
    double const* alpha,
    double const* A,
    double* B, lapack_int const* ldb );

#define LAPACK_stfsm LAPACK_GLOBAL(stfsm,STFSM)
void LAPACK_stfsm(
    char const* transr, char const* side, char const* uplo, char const* trans, char const* diag,
    lapack_int const* m, lapack_int const* n,
    float const* alpha,
    float const* A,
    float* B, lapack_int const* ldb );

#define LAPACK_ztfsm LAPACK_GLOBAL(ztfsm,ZTFSM)
void LAPACK_ztfsm(
    char const* transr, char const* side, char const* uplo, char const* trans, char const* diag,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* alpha,
    lapack_complex_double const* A,
    lapack_complex_double* B, lapack_int const* ldb );

#define LAPACK_ctftri LAPACK_GLOBAL(ctftri,CTFTRI)
void LAPACK_ctftri(
    char const* transr, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_float* A,
    lapack_int* info );

#define LAPACK_dtftri LAPACK_GLOBAL(dtftri,DTFTRI)
void LAPACK_dtftri(
    char const* transr, char const* uplo, char const* diag,
    lapack_int const* n,
    double* A,
    lapack_int* info );

#define LAPACK_stftri LAPACK_GLOBAL(stftri,STFTRI)
void LAPACK_stftri(
    char const* transr, char const* uplo, char const* diag,
    lapack_int const* n,
    float* A,
    lapack_int* info );

#define LAPACK_ztftri LAPACK_GLOBAL(ztftri,ZTFTRI)
void LAPACK_ztftri(
    char const* transr, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_double* A,
    lapack_int* info );

#define LAPACK_ctfttp LAPACK_GLOBAL(ctfttp,CTFTTP)
void LAPACK_ctfttp(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* ARF,
    lapack_complex_float* AP,
    lapack_int* info );

#define LAPACK_dtfttp LAPACK_GLOBAL(dtfttp,DTFTTP)
void LAPACK_dtfttp(
    char const* transr, char const* uplo,
    lapack_int const* n,
    double const* ARF,
    double* AP,
    lapack_int* info );

#define LAPACK_stfttp LAPACK_GLOBAL(stfttp,STFTTP)
void LAPACK_stfttp(
    char const* transr, char const* uplo,
    lapack_int const* n,
    float const* ARF,
    float* AP,
    lapack_int* info );

#define LAPACK_ztfttp LAPACK_GLOBAL(ztfttp,ZTFTTP)
void LAPACK_ztfttp(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* ARF,
    lapack_complex_double* AP,
    lapack_int* info );

#define LAPACK_ctfttr LAPACK_GLOBAL(ctfttr,CTFTTR)
void LAPACK_ctfttr(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* ARF,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dtfttr LAPACK_GLOBAL(dtfttr,DTFTTR)
void LAPACK_dtfttr(
    char const* transr, char const* uplo,
    lapack_int const* n,
    double const* ARF,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_stfttr LAPACK_GLOBAL(stfttr,STFTTR)
void LAPACK_stfttr(
    char const* transr, char const* uplo,
    lapack_int const* n,
    float const* ARF,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_ztfttr LAPACK_GLOBAL(ztfttr,ZTFTTR)
void LAPACK_ztfttr(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* ARF,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_ctgevc LAPACK_GLOBAL(ctgevc,CTGEVC)
void LAPACK_ctgevc(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_float const* S, lapack_int const* lds,
    lapack_complex_float const* P, lapack_int const* ldp,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtgevc LAPACK_GLOBAL(dtgevc,DTGEVC)
void LAPACK_dtgevc(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    double const* S, lapack_int const* lds,
    double const* P, lapack_int const* ldp,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    double* work,
    lapack_int* info );

#define LAPACK_stgevc LAPACK_GLOBAL(stgevc,STGEVC)
void LAPACK_stgevc(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    float const* S, lapack_int const* lds,
    float const* P, lapack_int const* ldp,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    float* work,
    lapack_int* info );

#define LAPACK_ztgevc LAPACK_GLOBAL(ztgevc,ZTGEVC)
void LAPACK_ztgevc(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_double const* S, lapack_int const* lds,
    lapack_complex_double const* P, lapack_int const* ldp,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctgexc LAPACK_GLOBAL(ctgexc,CTGEXC)
void LAPACK_ctgexc(
    lapack_logical const* wantq, lapack_logical const* wantz, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* Z, lapack_int const* ldz, lapack_int const* ifst, lapack_int* ilst,
    lapack_int* info );

#define LAPACK_dtgexc LAPACK_GLOBAL(dtgexc,DTGEXC)
void LAPACK_dtgexc(
    lapack_logical const* wantq, lapack_logical const* wantz, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* Q, lapack_int const* ldq,
    double* Z, lapack_int const* ldz, lapack_int* ifst, lapack_int* ilst,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_stgexc LAPACK_GLOBAL(stgexc,STGEXC)
void LAPACK_stgexc(
    lapack_logical const* wantq, lapack_logical const* wantz, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* Q, lapack_int const* ldq,
    float* Z, lapack_int const* ldz, lapack_int* ifst, lapack_int* ilst,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ztgexc LAPACK_GLOBAL(ztgexc,ZTGEXC)
void LAPACK_ztgexc(
    lapack_logical const* wantq, lapack_logical const* wantz, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* Z, lapack_int const* ldz, lapack_int const* ifst, lapack_int* ilst,
    lapack_int* info );

#define LAPACK_ctgsen LAPACK_GLOBAL(ctgsen,CTGSEN)
void LAPACK_ctgsen(
    lapack_int const* ijob, lapack_logical const* wantq, lapack_logical const* wantz, lapack_logical const* select, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* alpha,
    lapack_complex_float* beta,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* Z, lapack_int const* ldz, lapack_int* m,
    float* pl,
    float* pr,
    float* DIF,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_dtgsen LAPACK_GLOBAL(dtgsen,DTGSEN)
void LAPACK_dtgsen(
    lapack_int const* ijob, lapack_logical const* wantq, lapack_logical const* wantz, lapack_logical const* select, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* alphar,
    double* alphai,
    double* beta,
    double* Q, lapack_int const* ldq,
    double* Z, lapack_int const* ldz, lapack_int* m,
    double* pl,
    double* pr,
    double* DIF,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_stgsen LAPACK_GLOBAL(stgsen,STGSEN)
void LAPACK_stgsen(
    lapack_int const* ijob, lapack_logical const* wantq, lapack_logical const* wantz, lapack_logical const* select, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* alphar,
    float* alphai,
    float* beta,
    float* Q, lapack_int const* ldq,
    float* Z, lapack_int const* ldz, lapack_int* m,
    float* pl,
    float* pr,
    float* DIF,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ztgsen LAPACK_GLOBAL(ztgsen,ZTGSEN)
void LAPACK_ztgsen(
    lapack_int const* ijob, lapack_logical const* wantq, lapack_logical const* wantz, lapack_logical const* select, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* alpha,
    lapack_complex_double* beta,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* Z, lapack_int const* ldz, lapack_int* m,
    double* pl,
    double* pr,
    double* DIF,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ctgsja LAPACK_GLOBAL(ctgsja,CTGSJA)
void LAPACK_ctgsja(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    float const* tola,
    float const* tolb,
    float* alpha,
    float* beta,
    lapack_complex_float* U, lapack_int const* ldu,
    lapack_complex_float* V, lapack_int const* ldv,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* work, lapack_int* ncycle,
    lapack_int* info );

#define LAPACK_dtgsja LAPACK_GLOBAL(dtgsja,DTGSJA)
void LAPACK_dtgsja(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double const* tola,
    double const* tolb,
    double* alpha,
    double* beta,
    double* U, lapack_int const* ldu,
    double* V, lapack_int const* ldv,
    double* Q, lapack_int const* ldq,
    double* work, lapack_int* ncycle,
    lapack_int* info );

#define LAPACK_stgsja LAPACK_GLOBAL(stgsja,STGSJA)
void LAPACK_stgsja(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float const* tola,
    float const* tolb,
    float* alpha,
    float* beta,
    float* U, lapack_int const* ldu,
    float* V, lapack_int const* ldv,
    float* Q, lapack_int const* ldq,
    float* work, lapack_int* ncycle,
    lapack_int* info );

#define LAPACK_ztgsja LAPACK_GLOBAL(ztgsja,ZTGSJA)
void LAPACK_ztgsja(
    char const* jobu, char const* jobv, char const* jobq,
    lapack_int const* m, lapack_int const* p, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    double const* tola,
    double const* tolb,
    double* alpha,
    double* beta,
    lapack_complex_double* U, lapack_int const* ldu,
    lapack_complex_double* V, lapack_int const* ldv,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* work, lapack_int* ncycle,
    lapack_int* info );

#define LAPACK_ctgsna LAPACK_GLOBAL(ctgsna,CTGSNA)
void LAPACK_ctgsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float const* VL, lapack_int const* ldvl,
    lapack_complex_float const* VR, lapack_int const* ldvr,
    float* S,
    float* DIF, lapack_int const* mm, lapack_int* m,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dtgsna LAPACK_GLOBAL(dtgsna,DTGSNA)
void LAPACK_dtgsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double const* B, lapack_int const* ldb,
    double const* VL, lapack_int const* ldvl,
    double const* VR, lapack_int const* ldvr,
    double* S,
    double* DIF, lapack_int const* mm, lapack_int* m,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_stgsna LAPACK_GLOBAL(stgsna,STGSNA)
void LAPACK_stgsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float const* B, lapack_int const* ldb,
    float const* VL, lapack_int const* ldvl,
    float const* VR, lapack_int const* ldvr,
    float* S,
    float* DIF, lapack_int const* mm, lapack_int* m,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztgsna LAPACK_GLOBAL(ztgsna,ZTGSNA)
void LAPACK_ztgsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double const* VL, lapack_int const* ldvl,
    lapack_complex_double const* VR, lapack_int const* ldvr,
    double* S,
    double* DIF, lapack_int const* mm, lapack_int* m,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ctgsyl LAPACK_GLOBAL(ctgsyl,CTGSYL)
void LAPACK_ctgsyl(
    char const* trans,
    lapack_int const* ijob, lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float const* D, lapack_int const* ldd,
    lapack_complex_float const* E, lapack_int const* lde,
    lapack_complex_float* F, lapack_int const* ldf,
    float* dif,
    float* scale,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_dtgsyl LAPACK_GLOBAL(dtgsyl,DTGSYL)
void LAPACK_dtgsyl(
    char const* trans,
    lapack_int const* ijob, lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double const* B, lapack_int const* ldb,
    double* C, lapack_int const* ldc,
    double const* D, lapack_int const* ldd,
    double const* E, lapack_int const* lde,
    double* F, lapack_int const* ldf,
    double* dif,
    double* scale,
    double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_stgsyl LAPACK_GLOBAL(stgsyl,STGSYL)
void LAPACK_stgsyl(
    char const* trans,
    lapack_int const* ijob, lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float const* B, lapack_int const* ldb,
    float* C, lapack_int const* ldc,
    float const* D, lapack_int const* ldd,
    float const* E, lapack_int const* lde,
    float* F, lapack_int const* ldf,
    float* dif,
    float* scale,
    float* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztgsyl LAPACK_GLOBAL(ztgsyl,ZTGSYL)
void LAPACK_ztgsyl(
    char const* trans,
    lapack_int const* ijob, lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double const* D, lapack_int const* ldd,
    lapack_complex_double const* E, lapack_int const* lde,
    lapack_complex_double* F, lapack_int const* ldf,
    double* dif,
    double* scale,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ctpcon LAPACK_GLOBAL(ctpcon,CTPCON)
void LAPACK_ctpcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_float const* AP,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtpcon LAPACK_GLOBAL(dtpcon,DTPCON)
void LAPACK_dtpcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    double const* AP,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_stpcon LAPACK_GLOBAL(stpcon,STPCON)
void LAPACK_stpcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    float const* AP,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztpcon LAPACK_GLOBAL(ztpcon,ZTPCON)
void LAPACK_ztpcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_double const* AP,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctplqt LAPACK_GLOBAL(ctplqt,CTPLQT)
void LAPACK_ctplqt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* mb,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dtplqt LAPACK_GLOBAL(dtplqt,DTPLQT)
void LAPACK_dtplqt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* mb,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* T, lapack_int const* ldt,
    double* work,
    lapack_int* info );

#define LAPACK_stplqt LAPACK_GLOBAL(stplqt,STPLQT)
void LAPACK_stplqt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* mb,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* T, lapack_int const* ldt,
    float* work,
    lapack_int* info );

#define LAPACK_ztplqt LAPACK_GLOBAL(ztplqt,ZTPLQT)
void LAPACK_ztplqt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* mb,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_ctplqt2 LAPACK_GLOBAL(ctplqt2,CTPLQT2)
void LAPACK_ctplqt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_dtplqt2 LAPACK_GLOBAL(dtplqt2,DTPLQT2)
void LAPACK_dtplqt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_stplqt2 LAPACK_GLOBAL(stplqt2,STPLQT2)
void LAPACK_stplqt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_ztplqt2 LAPACK_GLOBAL(ztplqt2,ZTPLQT2)
void LAPACK_ztplqt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_ctpmlqt LAPACK_GLOBAL(ctpmlqt,CTPMLQT)
void LAPACK_ctpmlqt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* mb,
    lapack_complex_float const* V, lapack_int const* ldv,
    lapack_complex_float const* T, lapack_int const* ldt,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dtpmlqt LAPACK_GLOBAL(dtpmlqt,DTPMLQT)
void LAPACK_dtpmlqt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* mb,
    double const* V, lapack_int const* ldv,
    double const* T, lapack_int const* ldt,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* work,
    lapack_int* info );

#define LAPACK_stpmlqt LAPACK_GLOBAL(stpmlqt,STPMLQT)
void LAPACK_stpmlqt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* mb,
    float const* V, lapack_int const* ldv,
    float const* T, lapack_int const* ldt,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* work,
    lapack_int* info );

#define LAPACK_ztpmlqt LAPACK_GLOBAL(ztpmlqt,ZTPMLQT)
void LAPACK_ztpmlqt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* mb,
    lapack_complex_double const* V, lapack_int const* ldv,
    lapack_complex_double const* T, lapack_int const* ldt,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_ctpmqrt LAPACK_GLOBAL(ctpmqrt,CTPMQRT)
void LAPACK_ctpmqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* nb,
    lapack_complex_float const* V, lapack_int const* ldv,
    lapack_complex_float const* T, lapack_int const* ldt,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dtpmqrt LAPACK_GLOBAL(dtpmqrt,DTPMQRT)
void LAPACK_dtpmqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* nb,
    double const* V, lapack_int const* ldv,
    double const* T, lapack_int const* ldt,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* work,
    lapack_int* info );

#define LAPACK_stpmqrt LAPACK_GLOBAL(stpmqrt,STPMQRT)
void LAPACK_stpmqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* nb,
    float const* V, lapack_int const* ldv,
    float const* T, lapack_int const* ldt,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* work,
    lapack_int* info );

#define LAPACK_ztpmqrt LAPACK_GLOBAL(ztpmqrt,ZTPMQRT)
void LAPACK_ztpmqrt(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l, lapack_int const* nb,
    lapack_complex_double const* V, lapack_int const* ldv,
    lapack_complex_double const* T, lapack_int const* ldt,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_ctpqrt LAPACK_GLOBAL(ctpqrt,CTPQRT)
void LAPACK_ctpqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* nb,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_dtpqrt LAPACK_GLOBAL(dtpqrt,DTPQRT)
void LAPACK_dtpqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* nb,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* T, lapack_int const* ldt,
    double* work,
    lapack_int* info );

#define LAPACK_stpqrt LAPACK_GLOBAL(stpqrt,STPQRT)
void LAPACK_stpqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* nb,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* T, lapack_int const* ldt,
    float* work,
    lapack_int* info );

#define LAPACK_ztpqrt LAPACK_GLOBAL(ztpqrt,ZTPQRT)
void LAPACK_ztpqrt(
    lapack_int const* m, lapack_int const* n, lapack_int const* l, lapack_int const* nb,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_ctpqrt2 LAPACK_GLOBAL(ctpqrt2,CTPQRT2)
void LAPACK_ctpqrt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_dtpqrt2 LAPACK_GLOBAL(dtpqrt2,DTPQRT2)
void LAPACK_dtpqrt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_stpqrt2 LAPACK_GLOBAL(stpqrt2,STPQRT2)
void LAPACK_stpqrt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_ztpqrt2 LAPACK_GLOBAL(ztpqrt2,ZTPQRT2)
void LAPACK_ztpqrt2(
    lapack_int const* m, lapack_int const* n, lapack_int const* l,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_int* info );

#define LAPACK_ctprfb LAPACK_GLOBAL(ctprfb,CTPRFB)
void LAPACK_ctprfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    lapack_complex_float const* V, lapack_int const* ldv,
    lapack_complex_float const* T, lapack_int const* ldt,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_complex_float* work, lapack_int const* ldwork );

#define LAPACK_dtprfb LAPACK_GLOBAL(dtprfb,DTPRFB)
void LAPACK_dtprfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    double const* V, lapack_int const* ldv,
    double const* T, lapack_int const* ldt,
    double* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    double* work, lapack_int const* ldwork );

#define LAPACK_stprfb LAPACK_GLOBAL(stprfb,STPRFB)
void LAPACK_stprfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    float const* V, lapack_int const* ldv,
    float const* T, lapack_int const* ldt,
    float* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    float* work, lapack_int const* ldwork );

#define LAPACK_ztprfb LAPACK_GLOBAL(ztprfb,ZTPRFB)
void LAPACK_ztprfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    lapack_complex_double const* V, lapack_int const* ldv,
    lapack_complex_double const* T, lapack_int const* ldt,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_complex_double* work, lapack_int const* ldwork );

#define LAPACK_ctprfs LAPACK_GLOBAL(ctprfs,CTPRFS)
void LAPACK_ctprfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float const* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtprfs LAPACK_GLOBAL(dtprfs,DTPRFS)
void LAPACK_dtprfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    double const* AP,
    double const* B, lapack_int const* ldb,
    double const* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_stprfs LAPACK_GLOBAL(stprfs,STPRFS)
void LAPACK_stprfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    float const* AP,
    float const* B, lapack_int const* ldb,
    float const* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztprfs LAPACK_GLOBAL(ztprfs,ZTPRFS)
void LAPACK_ztprfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double const* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctptri LAPACK_GLOBAL(ctptri,CTPTRI)
void LAPACK_ctptri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_float* AP,
    lapack_int* info );

#define LAPACK_dtptri LAPACK_GLOBAL(dtptri,DTPTRI)
void LAPACK_dtptri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    double* AP,
    lapack_int* info );

#define LAPACK_stptri LAPACK_GLOBAL(stptri,STPTRI)
void LAPACK_stptri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    float* AP,
    lapack_int* info );

#define LAPACK_ztptri LAPACK_GLOBAL(ztptri,ZTPTRI)
void LAPACK_ztptri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_double* AP,
    lapack_int* info );

#define LAPACK_ctptrs LAPACK_GLOBAL(ctptrs,CTPTRS)
void LAPACK_ctptrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* AP,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dtptrs LAPACK_GLOBAL(dtptrs,DTPTRS)
void LAPACK_dtptrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    double const* AP,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_stptrs LAPACK_GLOBAL(stptrs,STPTRS)
void LAPACK_stptrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    float const* AP,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ztptrs LAPACK_GLOBAL(ztptrs,ZTPTRS)
void LAPACK_ztptrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* AP,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ctpttf LAPACK_GLOBAL(ctpttf,CTPTTF)
void LAPACK_ctpttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP,
    lapack_complex_float* ARF,
    lapack_int* info );

#define LAPACK_dtpttf LAPACK_GLOBAL(dtpttf,DTPTTF)
void LAPACK_dtpttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    double const* AP,
    double* ARF,
    lapack_int* info );

#define LAPACK_stpttf LAPACK_GLOBAL(stpttf,STPTTF)
void LAPACK_stpttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    float const* AP,
    float* ARF,
    lapack_int* info );

#define LAPACK_ztpttf LAPACK_GLOBAL(ztpttf,ZTPTTF)
void LAPACK_ztpttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP,
    lapack_complex_double* ARF,
    lapack_int* info );

#define LAPACK_ctpttr LAPACK_GLOBAL(ctpttr,CTPTTR)
void LAPACK_ctpttr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dtpttr LAPACK_GLOBAL(dtpttr,DTPTTR)
void LAPACK_dtpttr(
    char const* uplo,
    lapack_int const* n,
    double const* AP,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_stpttr LAPACK_GLOBAL(stpttr,STPTTR)
void LAPACK_stpttr(
    char const* uplo,
    lapack_int const* n,
    float const* AP,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_ztpttr LAPACK_GLOBAL(ztpttr,ZTPTTR)
void LAPACK_ztpttr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_ctrcon LAPACK_GLOBAL(ctrcon,CTRCON)
void LAPACK_ctrcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    float* rcond,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtrcon LAPACK_GLOBAL(dtrcon,DTRCON)
void LAPACK_dtrcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* rcond,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_strcon LAPACK_GLOBAL(strcon,STRCON)
void LAPACK_strcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* rcond,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztrcon LAPACK_GLOBAL(ztrcon,ZTRCON)
void LAPACK_ztrcon(
    char const* norm, char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    double* rcond,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctrevc LAPACK_GLOBAL(ctrevc,CTREVC)
void LAPACK_ctrevc(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtrevc LAPACK_GLOBAL(dtrevc,DTREVC)
void LAPACK_dtrevc(
    char const* side, char const* howmny,
    lapack_logical* select,
    lapack_int const* n,
    double const* T, lapack_int const* ldt,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    double* work,
    lapack_int* info );

#define LAPACK_strevc LAPACK_GLOBAL(strevc,STREVC)
void LAPACK_strevc(
    char const* side, char const* howmny,
    lapack_logical* select,
    lapack_int const* n,
    float const* T, lapack_int const* ldt,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    float* work,
    lapack_int* info );

#define LAPACK_ztrevc LAPACK_GLOBAL(ztrevc,ZTREVC)
void LAPACK_ztrevc(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctrevc3 LAPACK_GLOBAL(ctrevc3,CTREVC3)
void LAPACK_ctrevc3(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* VL, lapack_int const* ldvl,
    lapack_complex_float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_dtrevc3 LAPACK_GLOBAL(dtrevc3,DTREVC3)
void LAPACK_dtrevc3(
    char const* side, char const* howmny,
    lapack_logical* select,
    lapack_int const* n,
    double const* T, lapack_int const* ldt,
    double* VL, lapack_int const* ldvl,
    double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_strevc3 LAPACK_GLOBAL(strevc3,STREVC3)
void LAPACK_strevc3(
    char const* side, char const* howmny,
    lapack_logical* select,
    lapack_int const* n,
    float const* T, lapack_int const* ldt,
    float* VL, lapack_int const* ldvl,
    float* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ztrevc3 LAPACK_GLOBAL(ztrevc3,ZTREVC3)
void LAPACK_ztrevc3(
    char const* side, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* VL, lapack_int const* ldvl,
    lapack_complex_double* VR, lapack_int const* ldvr, lapack_int const* mm, lapack_int* m,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* info );

#define LAPACK_ctrexc LAPACK_GLOBAL(ctrexc,CTREXC)
void LAPACK_ctrexc(
    char const* compq,
    lapack_int const* n,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* Q, lapack_int const* ldq, lapack_int const* ifst, lapack_int const* ilst,
    lapack_int* info );

#define LAPACK_dtrexc LAPACK_GLOBAL(dtrexc,DTREXC)
void LAPACK_dtrexc(
    char const* compq,
    lapack_int const* n,
    double* T, lapack_int const* ldt,
    double* Q, lapack_int const* ldq, lapack_int* ifst, lapack_int* ilst,
    double* work,
    lapack_int* info );

#define LAPACK_strexc LAPACK_GLOBAL(strexc,STREXC)
void LAPACK_strexc(
    char const* compq,
    lapack_int const* n,
    float* T, lapack_int const* ldt,
    float* Q, lapack_int const* ldq, lapack_int* ifst, lapack_int* ilst,
    float* work,
    lapack_int* info );

#define LAPACK_ztrexc LAPACK_GLOBAL(ztrexc,ZTREXC)
void LAPACK_ztrexc(
    char const* compq,
    lapack_int const* n,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* Q, lapack_int const* ldq, lapack_int const* ifst, lapack_int const* ilst,
    lapack_int* info );

#define LAPACK_ctrrfs LAPACK_GLOBAL(ctrrfs,CTRRFS)
void LAPACK_ctrrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float const* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    lapack_complex_float* work,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtrrfs LAPACK_GLOBAL(dtrrfs,DTRRFS)
void LAPACK_dtrrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double const* B, lapack_int const* ldb,
    double const* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    double* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_strrfs LAPACK_GLOBAL(strrfs,STRRFS)
void LAPACK_strrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float const* B, lapack_int const* ldb,
    float const* X, lapack_int const* ldx,
    float* ferr,
    float* berr,
    float* work,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztrrfs LAPACK_GLOBAL(ztrrfs,ZTRRFS)
void LAPACK_ztrrfs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double const* X, lapack_int const* ldx,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctrsen LAPACK_GLOBAL(ctrsen,CTRSEN)
void LAPACK_ctrsen(
    char const* job, char const* compq,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_float* T, lapack_int const* ldt,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* W, lapack_int* m,
    float* s,
    float* sep,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dtrsen LAPACK_GLOBAL(dtrsen,DTRSEN)
void LAPACK_dtrsen(
    char const* job, char const* compq,
    lapack_logical const* select,
    lapack_int const* n,
    double* T, lapack_int const* ldt,
    double* Q, lapack_int const* ldq,
    double* WR,
    double* WI, lapack_int* m,
    double* s,
    double* sep,
    double* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_strsen LAPACK_GLOBAL(strsen,STRSEN)
void LAPACK_strsen(
    char const* job, char const* compq,
    lapack_logical const* select,
    lapack_int const* n,
    float* T, lapack_int const* ldt,
    float* Q, lapack_int const* ldq,
    float* WR,
    float* WI, lapack_int* m,
    float* s,
    float* sep,
    float* work, lapack_int const* lwork,
    lapack_int* iwork, lapack_int const* liwork,
    lapack_int* info );

#define LAPACK_ztrsen LAPACK_GLOBAL(ztrsen,ZTRSEN)
void LAPACK_ztrsen(
    char const* job, char const* compq,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_double* T, lapack_int const* ldt,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* W, lapack_int* m,
    double* s,
    double* sep,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ctrsna LAPACK_GLOBAL(ctrsna,CTRSNA)
void LAPACK_ctrsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_float const* T, lapack_int const* ldt,
    lapack_complex_float const* VL, lapack_int const* ldvl,
    lapack_complex_float const* VR, lapack_int const* ldvr,
    float* S,
    float* SEP, lapack_int const* mm, lapack_int* m,
    lapack_complex_float* work, lapack_int const* ldwork,
    float* rwork,
    lapack_int* info );

#define LAPACK_dtrsna LAPACK_GLOBAL(dtrsna,DTRSNA)
void LAPACK_dtrsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    double const* T, lapack_int const* ldt,
    double const* VL, lapack_int const* ldvl,
    double const* VR, lapack_int const* ldvr,
    double* S,
    double* SEP, lapack_int const* mm, lapack_int* m,
    double* work, lapack_int const* ldwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_strsna LAPACK_GLOBAL(strsna,STRSNA)
void LAPACK_strsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    float const* T, lapack_int const* ldt,
    float const* VL, lapack_int const* ldvl,
    float const* VR, lapack_int const* ldvr,
    float* S,
    float* SEP, lapack_int const* mm, lapack_int* m,
    float* work, lapack_int const* ldwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_ztrsna LAPACK_GLOBAL(ztrsna,ZTRSNA)
void LAPACK_ztrsna(
    char const* job, char const* howmny,
    lapack_logical const* select,
    lapack_int const* n,
    lapack_complex_double const* T, lapack_int const* ldt,
    lapack_complex_double const* VL, lapack_int const* ldvl,
    lapack_complex_double const* VR, lapack_int const* ldvr,
    double* S,
    double* SEP, lapack_int const* mm, lapack_int* m,
    lapack_complex_double* work, lapack_int const* ldwork,
    double* rwork,
    lapack_int* info );

#define LAPACK_ctrsyl LAPACK_GLOBAL(ctrsyl,CTRSYL)
void LAPACK_ctrsyl(
    char const* trana, char const* tranb,
    lapack_int const* isgn, lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* B, lapack_int const* ldb,
    lapack_complex_float* C, lapack_int const* ldc,
    float* scale,
    lapack_int* info );

#define LAPACK_dtrsyl LAPACK_GLOBAL(dtrsyl,DTRSYL)
void LAPACK_dtrsyl(
    char const* trana, char const* tranb,
    lapack_int const* isgn, lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double const* B, lapack_int const* ldb,
    double* C, lapack_int const* ldc,
    double* scale,
    lapack_int* info );

#define LAPACK_strsyl LAPACK_GLOBAL(strsyl,STRSYL)
void LAPACK_strsyl(
    char const* trana, char const* tranb,
    lapack_int const* isgn, lapack_int const* m, lapack_int const* n,
    float const* A, lapack_int const* lda,
    float const* B, lapack_int const* ldb,
    float* C, lapack_int const* ldc,
    float* scale,
    lapack_int* info );

#define LAPACK_ztrsyl LAPACK_GLOBAL(ztrsyl,ZTRSYL)
void LAPACK_ztrsyl(
    char const* trana, char const* tranb,
    lapack_int const* isgn, lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* C, lapack_int const* ldc,
    double* scale,
    lapack_int* info );

#define LAPACK_ctrtri LAPACK_GLOBAL(ctrtri,CTRTRI)
void LAPACK_ctrtri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_dtrtri LAPACK_GLOBAL(dtrtri,DTRTRI)
void LAPACK_dtrtri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_strtri LAPACK_GLOBAL(strtri,STRTRI)
void LAPACK_strtri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    float* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_ztrtri LAPACK_GLOBAL(ztrtri,ZTRTRI)
void LAPACK_ztrtri(
    char const* uplo, char const* diag,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_int* info );

#define LAPACK_ctrtrs LAPACK_GLOBAL(ctrtrs,CTRTRS)
void LAPACK_ctrtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_dtrtrs LAPACK_GLOBAL(dtrtrs,DTRTRS)
void LAPACK_dtrtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    double const* A, lapack_int const* lda,
    double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_strtrs LAPACK_GLOBAL(strtrs,STRTRS)
void LAPACK_strtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    float const* A, lapack_int const* lda,
    float* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ztrtrs LAPACK_GLOBAL(ztrtrs,ZTRTRS)
void LAPACK_ztrtrs(
    char const* uplo, char const* trans, char const* diag,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_ctrttf LAPACK_GLOBAL(ctrttf,CTRTTF)
void LAPACK_ctrttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* ARF,
    lapack_int* info );

#define LAPACK_dtrttf LAPACK_GLOBAL(dtrttf,DTRTTF)
void LAPACK_dtrttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* ARF,
    lapack_int* info );

#define LAPACK_strttf LAPACK_GLOBAL(strttf,STRTTF)
void LAPACK_strttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* ARF,
    lapack_int* info );

#define LAPACK_ztrttf LAPACK_GLOBAL(ztrttf,ZTRTTF)
void LAPACK_ztrttf(
    char const* transr, char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* ARF,
    lapack_int* info );

#define LAPACK_ctrttp LAPACK_GLOBAL(ctrttp,CTRTTP)
void LAPACK_ctrttp(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float* AP,
    lapack_int* info );

#define LAPACK_dtrttp LAPACK_GLOBAL(dtrttp,DTRTTP)
void LAPACK_dtrttp(
    char const* uplo,
    lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* AP,
    lapack_int* info );

#define LAPACK_strttp LAPACK_GLOBAL(strttp,STRTTP)
void LAPACK_strttp(
    char const* uplo,
    lapack_int const* n,
    float const* A, lapack_int const* lda,
    float* AP,
    lapack_int* info );

#define LAPACK_ztrttp LAPACK_GLOBAL(ztrttp,ZTRTTP)
void LAPACK_ztrttp(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double* AP,
    lapack_int* info );

#define LAPACK_ctzrzf LAPACK_GLOBAL(ctzrzf,CTZRZF)
void LAPACK_ctzrzf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_dtzrzf LAPACK_GLOBAL(dtzrzf,DTZRZF)
void LAPACK_dtzrzf(
    lapack_int const* m, lapack_int const* n,
    double* A, lapack_int const* lda,
    double* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_stzrzf LAPACK_GLOBAL(stzrzf,STZRZF)
void LAPACK_stzrzf(
    lapack_int const* m, lapack_int const* n,
    float* A, lapack_int const* lda,
    float* tau,
    float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_ztzrzf LAPACK_GLOBAL(ztzrzf,ZTZRZF)
void LAPACK_ztzrzf(
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunbdb LAPACK_GLOBAL(cunbdb,CUNBDB)
void LAPACK_cunbdb(
    char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    lapack_complex_float* X11, lapack_int const* ldx11,
    lapack_complex_float* X12, lapack_int const* ldx12,
    lapack_complex_float* X21, lapack_int const* ldx21,
    lapack_complex_float* X22, lapack_int const* ldx22,
    float* theta,
    float* phi,
    lapack_complex_float* TAUP1,
    lapack_complex_float* TAUP2,
    lapack_complex_float* TAUQ1,
    lapack_complex_float* TAUQ2,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunbdb LAPACK_GLOBAL(zunbdb,ZUNBDB)
void LAPACK_zunbdb(
    char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    lapack_complex_double* X11, lapack_int const* ldx11,
    lapack_complex_double* X12, lapack_int const* ldx12,
    lapack_complex_double* X21, lapack_int const* ldx21,
    lapack_complex_double* X22, lapack_int const* ldx22,
    double* theta,
    double* phi,
    lapack_complex_double* TAUP1,
    lapack_complex_double* TAUP2,
    lapack_complex_double* TAUQ1,
    lapack_complex_double* TAUQ2,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cuncsd LAPACK_GLOBAL(cuncsd,CUNCSD)
void LAPACK_cuncsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    lapack_complex_float* X11, lapack_int const* ldx11,
    lapack_complex_float* X12, lapack_int const* ldx12,
    lapack_complex_float* X21, lapack_int const* ldx21,
    lapack_complex_float* X22, lapack_int const* ldx22,
    float* theta,
    lapack_complex_float* U1, lapack_int const* ldu1,
    lapack_complex_float* U2, lapack_int const* ldu2,
    lapack_complex_float* V1T, lapack_int const* ldv1t,
    lapack_complex_float* V2T, lapack_int const* ldv2t,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zuncsd LAPACK_GLOBAL(zuncsd,ZUNCSD)
void LAPACK_zuncsd(
    char const* jobu1, char const* jobu2, char const* jobv1t, char const* jobv2t, char const* trans, char const* signs,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    lapack_complex_double* X11, lapack_int const* ldx11,
    lapack_complex_double* X12, lapack_int const* ldx12,
    lapack_complex_double* X21, lapack_int const* ldx21,
    lapack_complex_double* X22, lapack_int const* ldx22,
    double* theta,
    lapack_complex_double* U1, lapack_int const* ldu1,
    lapack_complex_double* U2, lapack_int const* ldu2,
    lapack_complex_double* V1T, lapack_int const* ldv1t,
    lapack_complex_double* V2T, lapack_int const* ldv2t,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cuncsd2by1 LAPACK_GLOBAL(cuncsd2by1,CUNCSD2BY1)
void LAPACK_cuncsd2by1(
    char const* jobu1, char const* jobu2, char const* jobv1t,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    lapack_complex_float* X11, lapack_int const* ldx11,
    lapack_complex_float* X21, lapack_int const* ldx21,
    float* theta,
    lapack_complex_float* U1, lapack_int const* ldu1,
    lapack_complex_float* U2, lapack_int const* ldu2,
    lapack_complex_float* V1T, lapack_int const* ldv1t,
    lapack_complex_float* work, lapack_int const* lwork,
    float* rwork, lapack_int const* lrwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_zuncsd2by1 LAPACK_GLOBAL(zuncsd2by1,ZUNCSD2BY1)
void LAPACK_zuncsd2by1(
    char const* jobu1, char const* jobu2, char const* jobv1t,
    lapack_int const* m, lapack_int const* p, lapack_int const* q,
    lapack_complex_double* X11, lapack_int const* ldx11,
    lapack_complex_double* X21, lapack_int const* ldx21,
    double* theta,
    lapack_complex_double* U1, lapack_int const* ldu1,
    lapack_complex_double* U2, lapack_int const* ldu2,
    lapack_complex_double* V1T, lapack_int const* ldv1t,
    lapack_complex_double* work, lapack_int const* lwork,
    double* rwork, lapack_int const* lrwork,
    lapack_int* iwork,
    lapack_int* info );

#define LAPACK_cungbr LAPACK_GLOBAL(cungbr,CUNGBR)
void LAPACK_cungbr(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zungbr LAPACK_GLOBAL(zungbr,ZUNGBR)
void LAPACK_zungbr(
    char const* vect,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunghr LAPACK_GLOBAL(cunghr,CUNGHR)
void LAPACK_cunghr(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunghr LAPACK_GLOBAL(zunghr,ZUNGHR)
void LAPACK_zunghr(
    lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunglq LAPACK_GLOBAL(cunglq,CUNGLQ)
void LAPACK_cunglq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunglq LAPACK_GLOBAL(zunglq,ZUNGLQ)
void LAPACK_zunglq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cungql LAPACK_GLOBAL(cungql,CUNGQL)
void LAPACK_cungql(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zungql LAPACK_GLOBAL(zungql,ZUNGQL)
void LAPACK_zungql(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cungqr LAPACK_GLOBAL(cungqr,CUNGQR)
void LAPACK_cungqr(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zungqr LAPACK_GLOBAL(zungqr,ZUNGQR)
void LAPACK_zungqr(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cungrq LAPACK_GLOBAL(cungrq,CUNGRQ)
void LAPACK_cungrq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zungrq LAPACK_GLOBAL(zungrq,ZUNGRQ)
void LAPACK_zungrq(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cungtr LAPACK_GLOBAL(cungtr,CUNGTR)
void LAPACK_cungtr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zungtr LAPACK_GLOBAL(zungtr,ZUNGTR)
void LAPACK_zungtr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmbr LAPACK_GLOBAL(cunmbr,CUNMBR)
void LAPACK_cunmbr(
    char const* vect, char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmbr LAPACK_GLOBAL(zunmbr,ZUNMBR)
void LAPACK_zunmbr(
    char const* vect, char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmhr LAPACK_GLOBAL(cunmhr,CUNMHR)
void LAPACK_cunmhr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmhr LAPACK_GLOBAL(zunmhr,ZUNMHR)
void LAPACK_zunmhr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* ilo, lapack_int const* ihi,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmlq LAPACK_GLOBAL(cunmlq,CUNMLQ)
void LAPACK_cunmlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmlq LAPACK_GLOBAL(zunmlq,ZUNMLQ)
void LAPACK_zunmlq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmql LAPACK_GLOBAL(cunmql,CUNMQL)
void LAPACK_cunmql(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmql LAPACK_GLOBAL(zunmql,ZUNMQL)
void LAPACK_zunmql(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmqr LAPACK_GLOBAL(cunmqr,CUNMQR)
void LAPACK_cunmqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmqr LAPACK_GLOBAL(zunmqr,ZUNMQR)
void LAPACK_zunmqr(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmrq LAPACK_GLOBAL(cunmrq,CUNMRQ)
void LAPACK_cunmrq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmrq LAPACK_GLOBAL(zunmrq,ZUNMRQ)
void LAPACK_zunmrq(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmrz LAPACK_GLOBAL(cunmrz,CUNMRZ)
void LAPACK_cunmrz(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmrz LAPACK_GLOBAL(zunmrz,ZUNMRZ)
void LAPACK_zunmrz(
    char const* side, char const* trans,
    lapack_int const* m, lapack_int const* n, lapack_int const* k, lapack_int const* l,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cunmtr LAPACK_GLOBAL(cunmtr,CUNMTR)
void LAPACK_cunmtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* A, lapack_int const* lda,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_zunmtr LAPACK_GLOBAL(zunmtr,ZUNMTR)
void LAPACK_zunmtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* A, lapack_int const* lda,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work, lapack_int const* lwork,
    lapack_int* info );

#define LAPACK_cupgtr LAPACK_GLOBAL(cupgtr,CUPGTR)
void LAPACK_cupgtr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_float const* AP,
    lapack_complex_float const* tau,
    lapack_complex_float* Q, lapack_int const* ldq,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zupgtr LAPACK_GLOBAL(zupgtr,ZUPGTR)
void LAPACK_zupgtr(
    char const* uplo,
    lapack_int const* n,
    lapack_complex_double const* AP,
    lapack_complex_double const* tau,
    lapack_complex_double* Q, lapack_int const* ldq,
    lapack_complex_double* work,
    lapack_int* info );

#define LAPACK_cupmtr LAPACK_GLOBAL(cupmtr,CUPMTR)
void LAPACK_cupmtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_float const* AP,
    lapack_complex_float const* tau,
    lapack_complex_float* C, lapack_int const* ldc,
    lapack_complex_float* work,
    lapack_int* info );

#define LAPACK_zupmtr LAPACK_GLOBAL(zupmtr,ZUPMTR)
void LAPACK_zupmtr(
    char const* side, char const* uplo, char const* trans,
    lapack_int const* m, lapack_int const* n,
    lapack_complex_double const* AP,
    lapack_complex_double const* tau,
    lapack_complex_double* C, lapack_int const* ldc,
    lapack_complex_double* work,
    lapack_int* info );

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LAPACK_H */
