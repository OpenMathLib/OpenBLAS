/***************************************************************************
 * Copyright (c) 2025, The OpenBLAS Project
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. Neither the name of the OpenBLAS project nor the names of
 * its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * *****************************************************************************/

#if defined(BFLOAT16) && defined(BFLOAT16CONVERSION)
static float
bfloat16tof32 (bfloat16 value)
{
  blasint one = 1;
  float result;
  sbf16tos_(&one, &value, &one, &result, &one);
  return result;
}

#ifdef BGEMM
static bfloat16 f32tobfloat16(float value) {
  blasint one = 1;
  bfloat16 result;
  sbstobf16_(&one, &value, &one, &result, &one);
  return result;
}
#endif

#ifdef BGEMM
#define ALPHA bfloat16tof32(alpha)
#define BETA bfloat16tof32(beta)
#define BF16TOF32(x) (bfloat16tof32(x))
#define F32TOBF16(x) (f32tobfloat16(x))
#else
#define ALPHA alpha
#define BETA beta
#define BF16TOF32(x) (bfloat16tof32(x))
#define F32TOBF16(x) x
#endif
#else
#define ALPHA alpha
#define BETA beta
#define BF16TOF32(x) x
#define F32TOBF16(x) x
#endif
