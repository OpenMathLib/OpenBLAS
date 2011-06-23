#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#define  GAM     4096.e0
#define  GAMSQ   16777216.e0
#define  RGAMSQ  5.9604645e-8

#ifdef DOUBLE
#define ABS(x) fabs(x)
#else
#define ABS(x) fabsf(x)
#endif

#ifndef CBLAS

void NAME(FLOAT *dd1, FLOAT *dd2, FLOAT *dx1, FLOAT *DY1, FLOAT *dparam){

  FLOAT	dy1 = *DY1;

#else

void CNAME(FLOAT *dd1, FLOAT *dd2, FLOAT *dx1, FLOAT dy1, FLOAT *dparam){

#endif

    FLOAT du, dp1, dp2, dq2, dq1, dh11, dh21, dh12, dh22;
    int igo, flag;
    FLOAT dtemp;

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#else
  PRINT_DEBUG_CNAME;
#endif

    dh11 = ZERO;
    dh12 = ZERO;
    dh21 = ZERO;
    dh22 = ZERO;

    if (*dd1 < ZERO) goto L60;

    dp2 = *dd2 * dy1;

    if (dp2 == ZERO) {
      flag = -2;
      goto L260;
    }

    dp1 = *dd1 * *dx1;
    dq2 =  dp2 * dy1;
    dq1 =  dp1 * *dx1;

    if (! (ABS(dq1) > ABS(dq2))) goto L40;

    dh21 = -(dy1) / *dx1;
    dh12 = dp2 / dp1;

    du = ONE - dh12 * dh21;

    if (du <= ZERO) goto L60;

    flag = 0;
    *dd1 /= du;
    *dd2 /= du;
    *dx1 *= du;

    goto L100;

L40:
    if (dq2 < ZERO) goto L60;

    flag = 1;
    dh11  = dp1 / dp2;
    dh22  = *dx1 / dy1;
    du    = ONE + dh11 * dh22;
    dtemp = *dd2 / du;
    *dd2  = *dd1 / du;
    *dd1  = dtemp;
    *dx1  = dy1 * du;
    goto L100;

L60:
    flag = -1;
    dh11 = ZERO;
    dh12 = ZERO;
    dh21 = ZERO;
    dh22 = ZERO;

    *dd1 = ZERO;
    *dd2 = ZERO;
    *dx1 = ZERO;
    goto L220;


L70:
    if (flag < 0) goto L90;
 
    if (flag > 0) goto L80;
 
    dh11 = ONE;
    dh22 = ONE;
    flag = -1;
    goto L90;

L80:
    dh21 = -ONE;
    dh12 = ONE;
    flag = -1;

L90:
    switch (igo) {
	case 0: goto L120;
	case 1: goto L150;
	case 2: goto L180;
	case 3: goto L210;
    }

L100:
    if (!(*dd1 <= RGAMSQ)) goto L130;
    if (*dd1 == ZERO) goto L160;
    igo = 0;
    goto L70;

L120:
    *dd1 *= GAM * GAM;
    *dx1 /= GAM;
    dh11 /= GAM;
    dh12 /= GAM;
    goto L100;

L130:
    if (! (*dd1 >= GAMSQ)) {
	goto L160;
    }
    igo = 1;
    goto L70;

L150:
    *dd1 /= GAM * GAM;
    *dx1 *= GAM;
    dh11 *= GAM;
    dh12 *= GAM;
    goto L130;

L160:
    if (! (ABS(*dd2) <= RGAMSQ)) {
	goto L190;
    }
    if (*dd2 == ZERO) {
	goto L220;
    }
    igo = 2;
    goto L70;

L180:
/* Computing 2nd power */
    *dd2 *= GAM * GAM;
    dh21 /= GAM;
    dh22 /= GAM;
    goto L160;

L190:
    if (! (ABS(*dd2) >= GAMSQ)) {
	goto L220;
    }
    igo = 3;
    goto L70;

L210:
/* Computing 2nd power */
    *dd2 /= GAM * GAM;
    dh21 *= GAM;
    dh22 *= GAM;
    goto L190;

L220:
    if (flag < 0) {
	goto L250;
    } else if (flag == 0) {
	goto L230;
    } else {
	goto L240;
    }
L230:
    dparam[2] = dh21;
    dparam[3] = dh12;
    goto L260;
L240:
    dparam[2] = dh11;
    dparam[4] = dh22;
    goto L260;
L250:
    dparam[1] = dh11;
    dparam[2] = dh21;
    dparam[3] = dh12;
    dparam[4] = dh22;
L260:
    dparam[0] = (FLOAT) flag;
    return;
}


