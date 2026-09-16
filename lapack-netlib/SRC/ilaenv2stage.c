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

/* > \brief \b ILAENV2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILAENV2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv2
stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv2
stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv2
stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILAENV2STAGE( ISPEC, NAME, OPTS, N1, N2, N3, N4 ) */

/*       CHARACTER*( * )    NAME, OPTS */
/*       INTEGER            ISPEC, N1, N2, N3, N4 */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ILAENV2STAGE is called from the LAPACK routines to choose problem-dependent */
/* > parameters for the local environment.  See ISPEC for a description of */
/* > the parameters. */
/* > It sets problem and machine dependent parameters useful for *_2STAGE and */
/* > related subroutines. */
/* > */
/* > ILAENV2STAGE returns an INTEGER */
/* > if ILAENV2STAGE >= 0: ILAENV2STAGE returns the value of the parameter */
/* >                       specified by ISPEC */
/* > if ILAENV2STAGE < 0:  if ILAENV2STAGE = -k, the k-th argument had an */
/* >                       illegal value. */
/* > */
/* > This version provides a set of parameters which should give good, */
/* > but not optimal, performance on many of the currently available */
/* > computers for the 2-stage solvers. Users are encouraged to modify this */
/* > subroutine to set the tuning parameters for their particular machine using */
/* > the option and problem size information in the arguments. */
/* > */
/* > This routine will not function correctly if it is converted to all */
/* > lower case.  Converting it to all upper case is allowed. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is INTEGER */
/* >          Specifies the parameter to be returned as the value of */
/* >          ILAENV2STAGE. */
/* >          = 1: the optimal blocksize nb for the reduction to BAND */
/* > */
/* >          = 2: the optimal blocksize ib for the eigenvectors */
/* >               singular vectors update routine */
/* > */
/* >          = 3: The length of the array that store the Housholder */
/* >               representation for the second stage */
/* >               Band to Tridiagonal or Bidiagonal */
/* > */
/* >          = 4: The workspace needed for the routine in input. */
/* > */
/* >          = 5: For future release. */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* >          NAME is CHARACTER*(*) */
/* >          The name of the calling subroutine, in either upper case or */
/* >          lower case. */
/* > \endverbatim */
/* > */
/* > \param[in] OPTS */
/* > \verbatim */
/* >          OPTS is CHARACTER*(*) */
/* >          The character options to the subroutine NAME, concatenated */
/* >          into a single character string.  For example, UPLO = 'U', */
/* >          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/* >          be specified as OPTS = 'UTN'. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* >          N2 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N3 */
/* > \verbatim */
/* >          N3 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N4 */
/* > \verbatim */
/* >          N4 is INTEGER */
/* >          Problem dimensions for the subroutine NAME; these may not all */
/* >          be required. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \author Nick R. Papior */

/* > \date July 2017 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The following conventions have been used when calling ILAENV2STAGE */
/* > from the LAPACK routines: */
/* >  1)  OPTS is a concatenation of all of the character options to */
/* >      subroutine NAME, in the same order that they appear in the */
/* >      argument list for NAME, even if they are not used in determining */
/* >      the value of the parameter specified by ISPEC. */
/* >  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
/* >      that they appear in the argument list for NAME.  N1 is used */
/* >      first, N2 second, and so on, and unused problem dimensions are */
/* >      passed a value of -1. */
/* >  3)  The parameter value returned by ILAENV2STAGE is checked for validity in */
/* >      the calling subroutine. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
integer ilaenv2stage_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer iparam2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *);
    integer iispec;


/*  -- LAPACK auxiliary routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     July 2017 */


/*  ===================================================================== */

    switch (*ispec) {
	case 1:  goto L10;
	case 2:  goto L10;
	case 3:  goto L10;
	case 4:  goto L10;
	case 5:  goto L10;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L10:

/*     2stage eigenvalues and SVD or related subroutines. */

    iispec = *ispec + 16;
    ret_val = iparam2stage_(&iispec, name__, opts, n1, n2, n3, n4);
    return ret_val;

/*     End of ILAENV2STAGE */

} /* ilaenv2stage_ */

