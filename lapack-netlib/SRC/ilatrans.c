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




/* > \brief \b ILATRANS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILATRANS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilatran
s.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilatran
s.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilatran
s.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILATRANS( TRANS ) */

/*       CHARACTER          TRANS */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine translates from a character string specifying a */
/* > transposition operation to the relevant BLAST-specified integer */
/* > constant. */
/* > */
/* > ILATRANS returns an INTEGER.  If ILATRANS < 0, then the input is not */
/* > a character indicating a transposition operator.  Otherwise ILATRANS */
/* > returns the constant value corresponding to TRANS. */
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
integer ilatrans_(char *trans)
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

    if (lsame_(trans, "N")) {
	ret_val = 111;
    } else if (lsame_(trans, "T")) {
	ret_val = 112;
    } else if (lsame_(trans, "C")) {
	ret_val = 113;
    } else {
	ret_val = -1;
    }
    return ret_val;

/*     End of ILATRANS */

} /* ilatrans_ */

