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
typedef blasint logical;

typedef char logical1;
typedef char integer1;
typedef int ftnlen;

#define TRUE_ (1)
#define FALSE_ (0)
#define i_nint(x) ((integer)u_nint(*(x)))
#define u_nint(__x) ((__x)>=0 ? floor((__x) + .5) : -floor(.5 - (__x)))
#define f2cmax(a,b) ((a) >= (b) ? (a) : (b))
#define f2cmin(a,b) ((a) <= (b) ? (a) : (b))

#define s_cat(lpp, rpp, rnp, np, llp) { 	ftnlen i, nc, ll; char *f__rp, *lp; 	ll = (llp); lp = (lpp); 	for(i=0; i < (int)*(np); ++i) {         	nc = ll; 	        if((rnp)[i] < nc) nc = (rnp)[i]; 	        ll -= nc;         	f__rp = (rpp)[i]; 	        while(--nc >= 0) *lp++ = *(f__rp)++;         } 	while(--ll >= 0) *lp++ = ' '; }
#define s_cmp(a,b,c,d) ((integer)strncmp((a),(b),f2cmin((c),(d))))
#define s_copy(A,B,C,D) { int __i,__m; for (__i=0, __m=f2cmin((C),(D)); __i<__m && (B)[__i] != 0; ++__i) (A)[__i] = (B)[__i]; }
#define sig_die(s, kill) { exit(1); }
#define s_stop(s, n) {exit(0);}
#define myexit_() break;

/* > \brief \b IPARMQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download IPARMQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.  -f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.  f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.  f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK ) */

/*       INTEGER            IHI, ILO, ISPEC, LWORK, N */
/*       CHARACTER          NAME*( * ), OPTS*( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >      This program sets problem and machine dependent parameters */
/* >      useful for xHSEQR and related subroutines for eigenvalue */
/* >      problems. It is called whenever */
/* >      IPARMQ is called with 12 <= ISPEC <= 16 */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is INTEGER */
/* >              ISPEC specifies which tunable parameter IPARMQ should */
/* >              return. */
/* > */
/* >              ISPEC=12: (INMIN)  Matrices of order nmin or less */
/* >                        are sent directly to xLAHQR, the implicit */
/* >                        double shift QR algorithm.  NMIN must be */
/* >                        at least 11. */
/* > */
/* >              ISPEC=13: (INWIN)  Size of the deflation window. */
/* >                        This is best set greater than or equal to */
/* >                        the number of simultaneous shifts NS. */
/* >                        Larger matrices benefit from larger deflation */
/* >                        windows. */
/* > */
/* >              ISPEC=14: (INIBL) Determines when to stop nibbling and */
/* >                        invest in an (expensive) multi-shift QR sweep. */
/* >                        If the aggressive early deflation subroutine */
/* >                        finds LD converged eigenvalues from an order */
/* >                        NW deflation window and LD > (NW*NIBBLE)/100, */
/* >                        then the next QR sweep is skipped and early */
/* >                        deflation is applied immediately to the */
/* >                        remaining active diagonal block.  Setting */
/* >                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a */
/* >                        multi-shift QR sweep whenever early deflation */
/* >                        finds a converged eigenvalue.  Setting */
/* >                        IPARMQ(ISPEC=14) greater than or equal to 100 */
/* >                        prevents TTQRE from skipping a multi-shift */
/* >                        QR sweep. */
/* > */
/* >              ISPEC=15: (NSHFTS) The number of simultaneous shifts in */
/* >                        a multi-shift QR iteration. */
/* > */
/* >              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the */
/* >                        following meanings. */
/* >                        0:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are not */
/* >                            accumulated when updating the */
/* >                            far-from-diagonal matrix entries. */
/* >                        1:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are */
/* >                            accumulated, and matrix-matrix */
/* >                            multiplication is used to update the */
/* >                            far-from-diagonal matrix entries. */
/* >                        2:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are */
/* >                            accumulated, and 2-by-2 block structure */
/* >                            is exploited during matrix-matrix */
/* >                            multiplies. */
/* >                        (If xTRMM is slower than xGEMM, then */
/* >                        IPARMQ(ISPEC=16)=1 may be more efficient than */
/* >                        IPARMQ(ISPEC=16)=2 despite the greater level of */
/* >                        arithmetic work implied by the latter choice.) */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* >          NAME is CHARACTER string */
/* >               Name of the calling subroutine */
/* > \endverbatim */
/* > */
/* > \param[in] OPTS */
/* > \verbatim */
/* >          OPTS is CHARACTER string */
/* >               This is a concatenation of the string arguments to */
/* >               TTQRE. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >               N is the order of the Hessenberg matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >               It is assumed that H is already upper triangular */
/* >               in rows and columns 1:ILO-1 and IHI+1:N. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >               The amount of workspace available. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >       Little is known about how best to choose these parameters. */
/* >       It is possible to use different values of the parameters */
/* >       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR. */
/* > */
/* >       It is probably best to choose different parameters for */
/* >       different matrices and different parameters at different */
/* >       times during the iteration, but this has not been */
/* >       implemented --- yet. */
/* > */
/* > */
/* >       The best choices of most of the parameters depend */
/* >       in an ill-understood way on the relative execution */
/* >       rate of xLAQR3 and xLAQR5 and on the nature of each */
/* >       particular eigenvalue problem.  Experiment may be the */
/* >       only practical way to determine which choices are most */
/* >       effective. */
/* > */
/* >       Following is a list of default values supplied by IPARMQ. */
/* >       These defaults may be adjusted in order to attain better */
/* >       performance in any particular computational environment. */
/* > */
/* >       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point. */
/* >                        Default: 75. (Must be at least 11.) */
/* > */
/* >       IPARMQ(ISPEC=13) Recommended deflation window size. */
/* >                        This depends on ILO, IHI and NS, the */
/* >                        number of simultaneous shifts returned */
/* >                        by IPARMQ(ISPEC=15).  The default for */
/* >                        (IHI-ILO+1) <= 500 is NS.  The default */
/* >                        for (IHI-ILO+1) > 500 is 3*NS/2. */
/* > */
/* >       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14. */
/* > */
/* >       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS. */
/* >                        a multi-shift QR iteration. */
/* > */
/* >                        If IHI-ILO+1 is ... */
/* > */
/* >                        greater than      ...but less    ... the */
/* >                        or equal to ...      than        default is */
/* > */
/* >                                0               30       NS =   2+ */
/* >                               30               60       NS =   4+ */
/* >                               60              150       NS =  10 */
/* >                              150              590       NS =  ** */
/* >                              590             3000       NS =  64 */
/* >                             3000             6000       NS = 128 */
/* >                             6000             infinity   NS = 256 */
/* > */
/* >                    (+)  By default matrices of this order are */
/* >                         passed to the implicit double shift routine */
/* >                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These */
/* >                         values of NS are used only in case of a rare */
/* >                         xLAHQR failure. */
/* > */
/* >                    (**) The asterisks (**) indicate an ad-hoc */
/* >                         function increasing from 10 to 64. */
/* > */
/* >       IPARMQ(ISPEC=16) Select structured matrix multiply. */
/* >                        (See ISPEC=16 above for details.) */
/* >                        Default: 3. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */

integer iparmq_(integer *ispec, char *name__, char *opts, integer *n, integer 
       *ilo, integer *ihi, integer *lwork)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1;

    /* Local variables */
    integer i__, ic, nh, ns, iz;
    char subnam[7];
    integer name_len=0;

/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */


/*  ================================================================ */
    if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

/*        ==== Set the number simultaneous shifts ==== */

	nh = *ihi - *ilo + 1;
	ns = 2;
	if (nh >= 30) {
	    ns = 4;
	}
	if (nh >= 60) {
	    ns = 10;
	}
	if (nh >= 150) {
/* Computing MAX */
	    r__1 = log((real) nh) / log(2.f);
	    i__1 = 10, i__2 = nh / i_nint(&r__1);
	    ns = f2cmax(i__1,i__2);
	}
	if (nh >= 590) {
	    ns = 64;
	}
	if (nh >= 3000) {
	    ns = 128;
	}
	if (nh >= 6000) {
	    ns = 256;
	}
/* Computing MAX */
	i__1 = 2, i__2 = ns - ns % 2;
	ns = f2cmax(i__1,i__2);
    }

    if (*ispec == 12) {


/*        ===== Matrices of order smaller than NMIN get sent */
/*        .     to xLAHQR, the classic double shift algorithm. */
/*        .     This must be at least 11. ==== */

	ret_val = 75;

    } else if (*ispec == 14) {

/*        ==== INIBL: skip a multi-shift qr iteration and */
/*        .    whenever aggressive early deflation finds */
/*        .    at least (NIBBLE*(window size)/100) deflations. ==== */

	ret_val = 14;

    } else if (*ispec == 15) {

/*        ==== NSHFTS: The number of simultaneous shifts ===== */

	ret_val = ns;

    } else if (*ispec == 13) {

/*        ==== NW: deflation window size.  ==== */

	if (nh <= 500) {
	    ret_val = ns;
	} else {
	    ret_val = ns * 3 / 2;
	}

    } else if (*ispec == 16) {

/*        ==== IACC22: Whether to accumulate reflections */
/*        .     before updating the far-from-diagonal elements */
/*        .     and whether to use 2-by-2 block structure while */
/*        .     doing it.  A small amount of work could be saved */
/*        .     by making this choice dependent also upon the */
/*        .     NH=IHI-ILO+1. */


/*        Convert NAME to upper case if the first character is lower case. */

	ret_val = 0;
//	s_copy(subnam, name__, (ftnlen)6, name_len);
	strncpy(subnam,name__,6);
        subnam[6]='\0';
	ic = *(unsigned char *)subnam;
	iz = 'Z';
	if (iz == 90 || iz == 122) {

/*           ASCII character set */

	    if (ic >= 97 && ic <= 122) {
		*(unsigned char *)subnam = (char) (ic - 32);
		for (i__ = 2; i__ <= 6; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 97 && ic <= 122) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		    }
		}
	    }

	} else if (iz == 233 || iz == 169) {

/*           EBCDIC character set */

	    if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 
		    && ic <= 169)) {
		*(unsigned char *)subnam = (char) (ic + 64);
		for (i__ = 2; i__ <= 6; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || 
			    (ic >= 162 && ic <= 169)) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		    }
		}
	    }

	} else if (iz == 218 || iz == 250) {

/*           Prime machines:  ASCII+128 */

	    if (ic >= 225 && ic <= 250) {
		*(unsigned char *)subnam = (char) (ic - 32);
		for (i__ = 2; i__ <= 6; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 225 && ic <= 250) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		    }
		}
	    }
	}

	if (s_cmp(subnam + 1, "GGHRD", (ftnlen)5, (ftnlen)5) == 0 || s_cmp(
		subnam + 1, "GGHD3", (ftnlen)5, (ftnlen)5) == 0) {
	    ret_val = 1;
	    if (nh >= 14) {
		ret_val = 2;
	    }
	} else if (s_cmp(subnam + 3, "EXC", (ftnlen)3, (ftnlen)3) == 0) {
	    if (nh >= 14) {
		ret_val = 1;
	    }
	    if (nh >= 14) {
		ret_val = 2;
	    }
	} else if (s_cmp(subnam + 1, "HSEQR", (ftnlen)5, (ftnlen)5) == 0 || 
		s_cmp(subnam + 1, "LAQR", (ftnlen)4, (ftnlen)4) == 0) {
	    if (ns >= 14) {
		ret_val = 1;
	    }
	    if (ns >= 14) {
		ret_val = 2;
	    }
	}

    } else {
/*        ===== invalid value of ispec ===== */
	ret_val = -1;

    }

/*     ==== End of IPARMQ ==== */

    return ret_val;
} /* iparmq_ */

