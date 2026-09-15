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
typedef blasint logical;

typedef unsigned int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;

/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/




/* > \brief \b ILAUPLO */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILAUPLO + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilauplo
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilauplo
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilauplo
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILAUPLO( UPLO ) */

/*       CHARACTER          UPLO */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine translated from a character string specifying a */
/* > upper- or lower-triangular matrix to the relevant BLAST-specified */
/* > integer constant. */
/* > */
/* > ILAUPLO returns an INTEGER.  If ILAUPLO < 0, then the input is not */
/* > a character indicating an upper- or lower-triangular matrix. */
/* > Otherwise ILAUPLO returns the constant value corresponding to UPLO. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */


/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
integer ilauplo_(char *uplo)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern logical lsame_(char *, char *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */


/*  ===================================================================== */

    if (lsame_(uplo, "U")) {
	ret_val = 121;
    } else if (lsame_(uplo, "L")) {
	ret_val = 122;
    } else {
	ret_val = -1;
    }
    return ret_val;

/*     End of ILAUPLO */

} /* ilauplo_ */

