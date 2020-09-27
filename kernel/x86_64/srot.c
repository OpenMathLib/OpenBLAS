#include "common.h"

#if defined(SKYLAKEX)
#include "srot_microk_skylakex-2.c"
#elif defined(HASWELL)
#include "srot_microk_haswell-2.c"
#endif

#ifndef HAVE_SROT_KERNEL

static void srot_kernel(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT c, FLOAT s)
{
    BLASLONG i = 0;
    FLOAT f0, f1, f2, f3;
    FLOAT x0, x1, x2, x3;
    FLOAT g0, g1, g2, g3;
    FLOAT y0, y1, y2, y3;

    FLOAT* xp = x;
    FLOAT* yp = y;

    BLASLONG n1 = n & (~7);

    while (i < n1) {
        x0 = xp[0];
        y0 = yp[0];
        x1 = xp[1];
        y1 = yp[1];
        x2 = xp[2];
        y2 = yp[2];
        x3 = xp[3];
        y3 = yp[3];

        f0 = c*x0 + s*y0;
        g0 = c*y0 - s*x0;
        f1 = c*x1 + s*y1;
        g1 = c*y1 - s*x1;
        f2 = c*x2 + s*y2;
        g2 = c*y2 - s*x2;
        f3 = c*x3 + s*y3;
        g3 = c*y3 - s*x3;

        xp[0] = f0;
        yp[0] = g0;
        xp[1] = f1;
        yp[1] = g1;
        xp[2] = f2;
        yp[2] = g2;
        xp[3] = f3;
        yp[3] = g3;

        xp += 4;
        yp += 4;
        i += 4;
    }

    while (i < n) {
        FLOAT temp = c*x[i] + s*y[i];
        y[i] = c*y[i] - s*x[i];
        x[i] = temp;

        i++;
    }
}

#endif
static void rot_compute(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s)
{
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;

    FLOAT temp;
    
    if (n <= 0)
        return;
    if ((inc_x == 1) && (inc_y == 1)) {
            srot_kernel(n, x, y, c, s);
    }
    else {
        while (i < n) {
            temp = c * x[ix] + s * y[iy];
            y[iy] = c * y[iy] - s * x[ix];
            x[ix] = temp;

            ix += inc_x;
            iy += inc_y;
            i++;
        }
    }
    return;
}


#if defined(SMP)
static int rot_thread_function(blas_arg_t *args)
{

    rot_compute(args->m, 
            args->a, args->lda, 
            args->b, args->ldb, 
            ((float *)args->alpha)[0], 
            ((float *)args->alpha)[1]);
    return 0;
}

extern int blas_level1_thread(int mode, BLASLONG m, BLASLONG n, BLASLONG k, void *alpha, void *a, BLASLONG lda, void *b, BLASLONG ldb, void *c, BLASLONG ldc, int (*function)(), int nthreads);
#endif
int CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s)
{
#if defined(SMP)
    int nthreads;
    FLOAT alpha[2]={c, s};
    FLOAT dummy_c;
#endif

#if defined(SMP)
    if (inc_x == 0 || inc_y == 0 || n <= 100000) {
        nthreads = 1;
    }
    else {
        nthreads = num_cpu_avail(1);
    }

    if (nthreads == 1) {
        rot_compute(n, x, inc_x, y, inc_y, c, s);
    }
    else {
#if defined(DOUBLE)
	    int mode = BLAS_DOUBLE | BLAS_REAL | BLAS_PTHREAD;
#else
	    int mode = BLAS_SINGLE | BLAS_REAL | BLAS_PTHREAD;
#endif
	    blas_level1_thread(mode, n, 0, 0, alpha, x, inc_x, y, inc_y, &dummy_c, 0, (void *)rot_thread_function, nthreads);
    }
#else	
    rot_compute(n, x, inc_x, y, inc_y, c, s);
#endif
    return 0;
}
