/***************************************************************************
Copyright (c) 2020, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT* alpha, FLOAT * B, BLASLONG ldb, FLOAT * C, BLASLONG ldc)
{
	FLOAT real, imag;
	int i, j, l;
	for(i = 0; i < M; i++){
		for(j = 0; j < N; j++){
			real=0;
			imag=0;

			for(l = 0; l < K; l++){
#if defined(TT)
				real += (A[i*2*lda + 2*l]*B[l*2*ldb + 2*j]
					 -A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j + 1]);

				imag+=(A[i*2*lda + 2*l] * B[l*2*ldb + 2*j + 1]
				       + A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j]);

#elif defined(TC)
				real += (A[i*2*lda + 2*l]*B[l*2*ldb + 2*j]
					 +A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j + 1]);

				imag+=(-A[i*2*lda + 2*l] * B[l*2*ldb + 2*j + 1]
				       + A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j]);

#elif defined(CT)
				real += (A[i*2*lda + 2*l]*B[l*2*ldb + 2*j]
					 +A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j + 1]);

				imag+=(A[i*2*lda + 2*l] * B[l*2*ldb + 2*j + 1]
				       - A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j]);

#elif defined(CC)
				real += (A[i*2*lda + 2*l]*B[l*2*ldb + 2*j]
					 -A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j + 1]);

				imag+=(-A[i*2*lda + 2*l] * B[l*2*ldb + 2*j + 1]
				       - A[i*2*lda + 2*l + 1] * B[l*2*ldb + 2*j]);

#endif
			}

			C[j*2*ldc + 2*i] = alpha[0]*real - alpha[1]*imag;
			C[j*2*ldc+ 2*i + 1] = alpha[0]*imag + real*alpha[1];
		}
	}
	
	return 0;
}
