#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>
#ifdef complex
#undef complex
#endif
#ifdef I
#undef I
#endif

#if defined(_WIN64)
typedef long long BLASLONG;
typedef unsigned long long BLASULONG;
#else
typedef long BLASLONG;
typedef unsigned long BLASULONG;
#endif

#ifdef LAPACK_ILP64
typedef BLASLONG blasint;
#if defined(_WIN64)
#define blasabs(x) llabs(x)
#else
#define blasabs(x) labs(x)
#endif
#else
typedef int blasint;
#define blasabs(x) abs(x)
#endif

typedef blasint integer;

typedef unsigned int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef blasint logical;

typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

#ifdef _MSC_VER
static inline _Fcomplex Cf(complex *z) {_Fcomplex zz={z->r , z->i}; return zz;}
#else
static inline _Complex float Cf(complex *z) {return z->r + z->i*_Complex_I;}
#endif

#define c_abs(z) (cabsf(Cf(z)))

/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/




/* > \brief \b ICMAX1 finds the index of the first vector element of maximum absolute value. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ICMAX1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/icmax1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/icmax1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/icmax1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER          FUNCTION ICMAX1( N, CX, INCX ) */

/*       INTEGER            INCX, N */
/*       COMPLEX            CX( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ICMAX1 finds the index of the first vector element of maximum absolute value. */
/* > */
/* > Based on ICAMAX from Level 1 BLAS. */
/* > The change is to use the 'genuine' absolute value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of elements in the vector CX. */
/* > \endverbatim */
/* > */
/* > \param[in] CX */
/* > \verbatim */
/* >          CX is COMPLEX array, dimension (N) */
/* >          The vector CX. The ICMAX1 function returns the index of its first */
/* >          element of maximum absolute value. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The spacing between successive values of CX.  INCX >= 1. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date February 2014 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Nick Higham for use with CLACON. */

/*  ===================================================================== */
integer icmax1_(integer *n, complex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    real smax;
    integer i__, ix;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     February 2014 */


/*  ===================================================================== */


    /* Parameter adjustments */
    --cx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {

/*        code for increment equal to 1 */

	smax = c_abs(&cx[1]);
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    if (c_abs(&cx[i__]) > smax) {
		ret_val = i__;
		smax = c_abs(&cx[i__]);
	    }
	}
    } else {

/*        code for increment not equal to 1 */

	ix = 1;
	smax = c_abs(&cx[1]);
	ix += *incx;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    if (c_abs(&cx[ix]) > smax) {
		ret_val = i__;
		smax = c_abs(&cx[ix]);
	    }
	    ix += *incx;
	}
    }
    return ret_val;

/*     End of ICMAX1 */

} /* icmax1_ */

