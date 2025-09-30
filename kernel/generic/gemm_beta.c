/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright 2025 The OpenBLAS Project.						         */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include "common.h"

#if defined(BFLOAT16) && defined(BGEMM) && defined(BFLOAT16CONVERSION)
static float
bfloat16tof32 (bfloat16 f16)
{
  float result = 0;
  unsigned short* q = (unsigned short*)(&result);
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  q[0] = f16;
#else
  q[1] = f16;
#endif
  return result;
}
static bfloat16
f32tobfloat16(float f32)
{
    unsigned short* q = (unsigned short*)(&f32);
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return q[0];
#else
    return q[1];
#endif
}

#define BF16TOF32(x) (bfloat16tof32(x))
#define F32TOBF16(x) (f32tobfloat16(x))
#else
#define BF16TOF32(x) x
#define F32TOBF16(x) x
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT beta,
	  IFLOAT *dummy2, BLASLONG dummy3, IFLOAT *dummy4, BLASLONG dummy5,
	  FLOAT *c, BLASLONG ldc){

  BLASLONG i, j;
  BLASLONG chunk, remain;
  FLOAT *c_offset1, *c_offset;
  c_offset = c;
  chunk = m >> 3;
  remain = m & 7;
  if (beta == ZERO){
	  for(j=n; j>0; j--){
		c_offset1 = c_offset;
		c_offset += ldc;
		for(i=chunk; i>0; i--){
			*(c_offset1 + 0) = F32TOBF16(ZERO);
			*(c_offset1 + 1) = F32TOBF16(ZERO);
			*(c_offset1 + 2) = F32TOBF16(ZERO);
			*(c_offset1 + 3) = F32TOBF16(ZERO);
			*(c_offset1 + 4) = F32TOBF16(ZERO);
			*(c_offset1 + 5) = F32TOBF16(ZERO);
			*(c_offset1 + 6) = F32TOBF16(ZERO);
			*(c_offset1 + 7) = F32TOBF16(ZERO);
			c_offset1 += 8;
		}
		for(i=remain; i>0; i--){
			*c_offset1 = F32TOBF16(ZERO);
			c_offset1 ++;
		}
	  }
  } else {
	  for(j=n; j>0; j--){
		c_offset1 = c_offset;
		c_offset += ldc;
		for(i=chunk; i>0; i--){
			*(c_offset1 + 0) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[0]));
			*(c_offset1 + 1) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[1]));
			*(c_offset1 + 2) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[2]));
			*(c_offset1 + 3) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[3]));
			*(c_offset1 + 4) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[4]));
			*(c_offset1 + 5) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[5]));
			*(c_offset1 + 6) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[6]));
			*(c_offset1 + 7) = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[7]));
			c_offset1 += 8;
		}
		for(i=remain; i>0; i--){
			*c_offset1 = F32TOBF16(BF16TOF32(beta) * BF16TOF32(c_offset1[0]));
			c_offset1 ++;
		}
	  }
  }
  return 0;
};
