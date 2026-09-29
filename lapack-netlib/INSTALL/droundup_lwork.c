#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

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
typedef int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/*      INTEGER          LWORK */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DROUNDUP_LWORK deals with a subtle bug with returning LWORK as a Float. */
/* > This routine guarantees it is rounded up instead of down by */
/* > multiplying LWORK by 1+eps when it is necessary, where eps is the relative machine precision. */
/* > E.g., */
/* > */
/* >        float( 9007199254740993            ) == 9007199254740992 */
/* >        float( 9007199254740993 ) * (1.+eps) == 9007199254740994 */
/* > */
/* > \return DROUNDUP_LWORK */
/* > \verbatim */
/* >         DROUNDUP_LWORK >= LWORK. */
/* >         DROUNDUP_LWORK is guaranteed to have zero decimal part. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] LWORK Workspace size. */

/*  Authors: */
/*  ======== */

/* > \author Weslley Pereira, University of Colorado Denver, USA */

/* > \ingroup auxOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* >  This routine was inspired in the method `magma_zmake_lwork` from MAGMA. */
/* >  \see https://bitbucket.org/icl/magma/src/master/control/magma_zauxiliary.cpp */
/* > \endverbatim */

/*  ===================================================================== */
double droundup_lwork__(integer *lwork)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double epsilon=DBL_EPSILON;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/* ===================================================================== */
    ret_val = (double) (*lwork);

    if ((integer) ret_val < *lwork) {
/*         Force round up of LWORK */
	ret_val *= epsilon + 1.f;
    }

    return ret_val;

/*     End of DROUNDUP_LWORK */

} /* sroundup_lwork__ */

