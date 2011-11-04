#include "common.h"

//These are auto-tuning codes on Loongson-3A platform. 
//#define prefetch(x) __builtin_prefetch(x)
//#define prefetch(x) do {_mm_prefetch((char *)(x), _MM_HINT_T0);} while(0)
#define prefetch(x) __asm__ __volatile__("ld $0, %0"::"m"(x))
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define spec_loop_alpha1 do {Y[k] += A[jj + ii] * X[ii]; Y[k + 1] += A[jj + ii + 1] * X[ii]; Y[k + 1] += A[jj + ii] * X[ii + 1]; Y[k] -= A[jj + ii + 1] * X[ii + 1]; ii += 2;} while(0)
//#define spec_loop_alpha1 do {Y[ii] += A[jj + ii] * X[k] - A[jj + ii + 1] * X[k + 1]; Y[ii + 1] += A[jj + ii + 1] * X[k] + A[jj + ii] * X[k + 1]; ii += 2;} while(0)
#define spec_loop do {rTmp = A[jj + ii] * X[ii] - A[jj + ii + 1] * X[ii + 1]; iTmp = A[jj + ii] * X[ii + 1] + A[jj + ii + 1] * X[ii]; Y[k] += rTmp * rALPHA - iTmp * iALPHA; Y[k + 1] += rTmp * iALPHA + iTmp * rALPHA; ii += 2;} while(0)
#define norm_loop_alpha1 do {Y[k] += A[jj + ii] * X[iii] - A[jj + ii + 1] * X[iii + 1]; Y[k + 1] += A[jj + ii] * X[iii + 1] + A[jj + ii + 1] * X[iii]; ii += 2; iii += INCX * 2;} while(0)
#define norm_loop do {rTmp = A[jj + ii] * X[iii] - A[jj + ii + 1] * X[iii + 1]; iTmp = A[jj + ii] * X[iii + 1] + A[jj + ii + 1] * X[iii]; Y[k] += rTmp * rALPHA - iTmp * iALPHA; Y[k + 1] += rTmp * iALPHA + iTmp * rALPHA; ii += 2; iii += INCX * 2;} while(0)

int CNAME(BLASLONG M, BLASLONG N, BLASLONG UNUSED, FLOAT rALPHA, FLOAT iALPHA, FLOAT *A, BLASLONG LDA, FLOAT *X, BLASLONG INCX, FLOAT *Y, BLASLONG INCY, FLOAT *BUFFER) {

	if(!rALPHA && iALPHA)
		return 0;

//	if(INCX < 0)
//		INCX = -INCX;
//	if(INCY < 0)
//		INCY = -INCY;

	BLASLONG fahead = 30;
	BLASLONG spec_unroll = 2;
	BLASLONG tMQ = M - M % spec_unroll;
	BLASLONG j = 0, k = 0, jj=0;


	if(rALPHA == 1 && iALPHA == 0) {
		if(INCX == 1) {
			for(; likely(j < N); j++, k += INCY * 2, jj += LDA * 2) {
				BLASLONG i = 0, ii = 0;
				for(; likely(i < tMQ); i += spec_unroll) {
					prefetch(A[jj + ii + fahead]);
					prefetch(X[ii + fahead]);
					/*loop_mark*/ spec_loop_alpha1;
					/*loop_mark*/ spec_loop_alpha1;
				}
				for(; likely(i < M); i++) {
					spec_loop_alpha1;
				}
			}
		} else {
			for(; likely(j < N); j++, k += INCY * 2, jj += LDA * 2) {
				BLASLONG i = 0, ii = 0, iii = 0;
				for(; likely(i < tMQ); i += spec_unroll) {
					prefetch(A[jj + ii + fahead]);
					prefetch(X[iii + fahead]);
					/*loop_mark*/ norm_loop_alpha1;
					/*loop_mark*/ norm_loop_alpha1;
				}
				for(; likely(i < M); i++) {
					norm_loop_alpha1;
				}
			}
		}
	} else {
		FLOAT rTmp, iTmp;
		if(INCX == 1) {
			for(; likely(j < N); j++, k += INCY * 2, jj += LDA * 2) {
				BLASLONG i = 0, ii = 0;
				for(; likely(i < tMQ); i += spec_unroll) {
					prefetch(A[jj + ii + fahead]);
					prefetch(X[ii + fahead]);
					/*loop_mark*/ spec_loop;
					/*loop_mark*/ spec_loop;
				}
				for(; likely(i < M); i++) {
					spec_loop;
				}
			}
		} else {
			for(; likely(j < N); j++, k += INCY * 2, jj += LDA * 2) {
				BLASLONG i = 0, ii = 0, iii = 0;
				for(; likely(i < tMQ); i += spec_unroll) {
					prefetch(A[jj + ii + fahead]);
					prefetch(X[iii + fahead]);
					/*loop_mark*/ norm_loop;
					/*loop_mark*/ norm_loop;
				}
				for(; likely(i < M); i++) {
					norm_loop;
				}
			}
		}
	}
	return 0;
}
