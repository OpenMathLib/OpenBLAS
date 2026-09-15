#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
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

/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b IPARAM2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download IPARAM2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparam2
stage.F"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparam2
stage.F"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparam2
stage.F"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IPARAM2STAGE( ISPEC, NAME, OPTS, */
/*                                    NI, NBI, IBI, NXI ) */
/*       #if defined(_OPENMP) */
/*           use omp_lib */
/*       #endif */
/*       IMPLICIT NONE */

/*       CHARACTER*( * )    NAME, OPTS */
/*       INTEGER            ISPEC, NI, NBI, IBI, NXI */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >      This program sets problem and machine dependent parameters */
/* >      useful for xHETRD_2STAGE, xHETRD_HE2HB, xHETRD_HB2ST, */
/* >      xGEBRD_2STAGE, xGEBRD_GE2GB, xGEBRD_GB2BD */
/* >      and related subroutines for eigenvalue problems. */
/* >      It is called whenever ILAENV is called with 17 <= ISPEC <= 21. */
/* >      It is called whenever ILAENV2STAGE is called with 1 <= ISPEC <= 5 */
/* >      with a direct conversion ISPEC + 16. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is integer scalar */
/* >              ISPEC specifies which tunable parameter IPARAM2STAGE should */
/* >              return. */
/* > */
/* >              ISPEC=17: the optimal blocksize nb for the reduction to */
/* >                        BAND */
/* > */
/* >              ISPEC=18: the optimal blocksize ib for the eigenvectors */
/* >                        singular vectors update routine */
/* > */
/* >              ISPEC=19: The length of the array that store the Housholder */
/* >                        representation for the second stage */
/* >                        Band to Tridiagonal or Bidiagonal */
/* > */
/* >              ISPEC=20: The workspace needed for the routine in input. */
/* > */
/* >              ISPEC=21: For future release. */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* >          NAME is character string */
/* >               Name of the calling subroutine */
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
/* > \param[in] NI */
/* > \verbatim */
/* >          NI is INTEGER which is the size of the matrix */
/* > \endverbatim */
/* > */
/* > \param[in] NBI */
/* > \verbatim */
/* >          NBI is INTEGER which is the used in the reduciton, */
/* >          (e.g., the size of the band), needed to compute workspace */
/* >          and LHOUS2. */
/* > \endverbatim */
/* > */
/* > \param[in] IBI */
/* > \verbatim */
/* >          IBI is INTEGER which represent the IB of the reduciton, */
/* >          needed to compute workspace and LHOUS2. */
/* > \endverbatim */
/* > */
/* > \param[in] NXI */
/* > \verbatim */
/* >          NXI is INTEGER needed in the future release. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Implemented by Azzam Haidar. */
/* > */
/* >  All detail are available on technical report, SC11, SC13 papers. */
/* > */
/* >  Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* >  Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* >  using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* >  of 2011 International Conference for High Performance Computing, */
/* >  Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* >  Article 8 , 11 pages. */
/* >  http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* >  A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* >  An improved parallel singular value algorithm and its implementation */
/* >  for multicore hardware, In Proceedings of 2013 International Conference */
/* >  for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* >  Denver, Colorado, USA, 2013. */
/* >  Article 90, 12 pages. */
/* >  http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* >  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* >  A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* >  calculations based on fine-grained memory aware tasks. */
/* >  International Journal of High Performance Computing Applications. */
/* >  Volume 28 Issue 2, Pages 196-209, May 2014. */
/* >  http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
integer iparam2stage_(integer *ispec, char *name__, char *opts, integer *ni, 
	integer *nbi, integer *ibi, integer *nxi)
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3;

    /* Local variables */
    char algo[4], prec[1], stag[6], vect[1];
    integer nthreads, i__;
    logical cprec, rprec;
    integer lhous, lwork, factoptnb, ib, ic, kd, iz;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    char subnam[14];
    integer lqoptnb, qroptnb;


/*  -- LAPACK auxiliary routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */


/*  ================================================================ */

/*     Invalid value for ISPEC */

    if (*ispec < 17 || *ispec > 21) {
	ret_val = -1;
	return ret_val;
    }

/*     Get the number of threads */

    nthreads = 1;
/*      WRITE(*,*) 'IPARAM VOICI NTHREADS ISPEC ',NTHREADS, ISPEC */

    if (*ispec != 19) {

/*        Convert NAME to upper case if the first character is lower case. */

	ret_val = -1;

//	s_copy(subnam, name__, (ftnlen)12, name_len);
	strncpy(subnam,name__,13);
	subnam[13]='\0';
	{
	    int i;
	    for (i=0;i<13;i++) subnam[i]=toupper(subnam[i]);
	}	

#if 0

	ic = *(unsigned char *)subnam;
	iz = 'Z';
	if (iz == 90 || iz == 122) {

/*           ASCII character set */

	    if (ic >= 97 && ic <= 122) {
		*(unsigned char *)subnam = (char) (ic - 32);
		for (i__ = 2; i__ <= 12; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 97 && ic <= 122) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		    }
/* L100: */
		}
	    }

	} else if (iz == 233 || iz == 169) {

/*           EBCDIC character set */

	    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 
		    && ic <= 169) {
		*(unsigned char *)subnam = (char) (ic + 64);
		for (i__ = 2; i__ <= 12; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || 
			    ic >= 162 && ic <= 169) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		    }
/* L110: */
		}
	    }

	} else if (iz == 218 || iz == 250) {

/*           Prime machines:  ASCII+128 */

	    if (ic >= 225 && ic <= 250) {
		*(unsigned char *)subnam = (char) (ic - 32);
		for (i__ = 2; i__ <= 12; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 225 && ic <= 250) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		    }
/* L120: */
		}
	    }
	}
#endif

//fprintf(stderr,"iparam2stage, subnam gross #%s#\n",subnam);

//	*(unsigned char *)prec = *(unsigned char *)subnam;
	strncpy(prec,subnam,1);
	strncpy(algo, subnam+3,3);
	algo[3]='\0';
	strncpy(stag, subnam+7,5);
	stag[5]='\0';
//	s_copy(algo, subnam + 3, (ftnlen)3, (ftnlen)3);
//	s_copy(stag, subnam + 7, (ftnlen)5, (ftnlen)5);
	rprec = *(unsigned char *)prec == 'S' || *(unsigned char *)prec == 
		'D';
	cprec = *(unsigned char *)prec == 'C' || *(unsigned char *)prec == 
		'Z';

/*        Invalid value for PRECISION */
//fprintf(stderr," prec %s algo %s stag %s\n",prec,algo,stag);
	if (! (rprec || cprec)) {
	    ret_val = -1;
	    return ret_val;
	}
    }
/*      WRITE(*,*),'RPREC,CPREC ',RPREC,CPREC, */
/*     $           '   ALGO ',ALGO,'    STAGE ',STAG */


    if (*ispec == 17 || *ispec == 18) {
//fprintf(stderr,"iparam2stage spec 17/18");
/*     ISPEC = 17, 18:  block size KD, IB */
/*     Could be also dependent from N but for now it */
/*     depend only on sequential or parallel */

	if (nthreads > 4) {
	    if (cprec) {
		kd = 128;
		ib = 32;
	    } else {
		kd = 160;
		ib = 40;
	    }
	} else if (nthreads > 1) {
	    if (cprec) {
		kd = 64;
		ib = 32;
	    } else {
		kd = 64;
		ib = 32;
	    }
	} else {
	    if (cprec) {
		kd = 16;
		ib = 16;
	    } else {
		kd = 32;
		ib = 16;
	    }
	}
	if (*ispec == 17) {
	    ret_val = kd;
	}
	if (*ispec == 18) {
	    ret_val = ib;
	}

    } else if (*ispec == 19) {
//fprintf(stderr,"iparam2stage spec 19\n");
/*     ISPEC = 19: */
/*     LHOUS length of the Houselholder representation */
/*     matrix (V,T) of the second stage. should be >= 1. */

/*     Will add the VECT OPTION HERE next release */
	*(unsigned char *)vect = *(unsigned char *)opts;
	if (*(unsigned char *)vect == 'N') {
/* Computing MAX */
	    i__1 = 1, i__2 = *ni << 2;
	    lhous = f2cmax(i__1,i__2);
	} else {
/*           This is not correct, it need to call the ALGO and the stage2 */
/* Computing MAX */
	    i__1 = 1, i__2 = *ni << 2;
	    lhous = f2cmax(i__1,i__2) + *ibi;
	}
	if (lhous >= 0) {
	    ret_val = lhous;
	} else {
	    ret_val = -1;
	}

    } else if (*ispec == 20) {
//fprintf(stderr,"iparam2stage spec 20\n");
/*     ISPEC = 20: (21 for future use) */
/*     LWORK length of the workspace for */
/*     either or both stages for TRD and BRD. should be >= 1. */
/*     TRD: */
/*     TRD_stage 1: = LT + LW + LS1 + LS2 */
/*                  = LDT*KD + N*KD + N*MAX(KD,FACTOPTNB) + LDS2*KD */
/*                    where LDT=LDS2=KD */
/*                  = N*KD + N*f2cmax(KD,FACTOPTNB) + 2*KD*KD */
/*     TRD_stage 2: = (2NB+1)*N + KD*NTHREADS */
/*     TRD_both   : = f2cmax(stage1,stage2) + AB ( AB=(KD+1)*N ) */
/*                  = N*KD + N*f2cmax(KD+1,FACTOPTNB) */
/*                    + f2cmax(2*KD*KD, KD*NTHREADS) */
/*                    + (KD+1)*N */
	lwork = -1;
	char *subnam=malloc(7*sizeof(char));
			strncpy(subnam,prec,1);
			sprintf(subnam+1,"GEQRF");
//	*(unsigned char *)subnam = *(unsigned char *)prec;
//	s_copy(subnam + 1, "GEQRF", (ftnlen)5, (ftnlen)5);
	qroptnb = ilaenv_(&c__1, subnam, " ", ni, nbi, &c_n1, &c_n1, (ftnlen)
		12, (ftnlen)1);
	sprintf(subnam+1,"GELQF");
	s_copy(subnam + 1, "GELQF", (ftnlen)5, (ftnlen)5);
	lqoptnb = ilaenv_(&c__1, subnam, " ", nbi, ni, &c_n1, &c_n1, (ftnlen)
		12, (ftnlen)1);
/*        Could be QR or LQ for TRD and the f2cmax for BRD */
	factoptnb = f2cmax(qroptnb,lqoptnb);
	if (s_cmp(algo, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (s_cmp(stag, "2STAG", (ftnlen)5, (ftnlen)5) == 0) {
/* Computing MAX */
		i__1 = *nbi + 1;
/* Computing MAX */
		i__2 = (*nbi << 1) * *nbi, i__3 = *nbi * nthreads;
		lwork = *ni * *nbi + *ni * f2cmax(i__1,factoptnb) + f2cmax(i__2,
			i__3) + (*nbi + 1) * *ni;
	    } else if (s_cmp(stag, "HE2HB", (ftnlen)5, (ftnlen)5) == 0 || 
		    s_cmp(stag, "SY2SB", (ftnlen)5, (ftnlen)5) == 0) {
		lwork = *ni * *nbi + *ni * f2cmax(*nbi,factoptnb) + (*nbi << 1) *
			 *nbi;
	    } else if (s_cmp(stag, "HB2ST", (ftnlen)5, (ftnlen)5) == 0 || 
		    s_cmp(stag, "SB2ST", (ftnlen)5, (ftnlen)5) == 0) {
		lwork = ((*nbi << 1) + 1) * *ni + *nbi * nthreads;
	    }
	} else if (s_cmp(algo, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (s_cmp(stag, "2STAG", (ftnlen)5, (ftnlen)5) == 0) {
/* Computing MAX */
		i__1 = *nbi + 1;
/* Computing MAX */
		i__2 = (*nbi << 1) * *nbi, i__3 = *nbi * nthreads;
		lwork = (*ni << 1) * *nbi + *ni * f2cmax(i__1,factoptnb) + f2cmax(
			i__2,i__3) + (*nbi + 1) * *ni;
	    } else if (s_cmp(stag, "GE2GB", (ftnlen)5, (ftnlen)5) == 0) {
		lwork = *ni * *nbi + *ni * f2cmax(*nbi,factoptnb) + (*nbi << 1) *
			 *nbi;
	    } else if (s_cmp(stag, "GB2BD", (ftnlen)5, (ftnlen)5) == 0) {
		lwork = (*nbi * 3 + 1) * *ni + *nbi * nthreads;
	    }
	}
	lwork = f2cmax(1,lwork);
	if (lwork > 0) {
	    ret_val = lwork;
	} else {
	    ret_val = -1;
	}

    } else if (*ispec == 21) {
//fprintf(stderr,"iparam2stage spec 21\n");
/*     ISPEC = 21 for future use */
	ret_val = *nxi;
    }

/*     ==== End of IPARAM2STAGE ==== */

    return ret_val;
} /* iparam2stage_ */

