/***************************************************************************
Copyright (c) 2025 The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H
#include <stdbool.h>

#include "../common.h"

#ifdef IBFLOAT16
static float float16to32(bfloat16 value)
{
  blasint one = 1;
  float result;
  sbf16tos_(&one, &value, &one, &result, &one);
  return result;
}
#endif

#ifdef OBFLOAT16
static float truncate_float32_to_bfloat16(float value) {
  blasint one = 1;
  bfloat16 tmp;
  float result;
  sbstobf16_(&one, &value, &one, &tmp, &one);
  sbf16tos_(&one, &tmp, &one, &result, &one);
  return result;
}
#endif

static void *malloc_safe(size_t size) {
  if (size == 0)
    return malloc(1);
  else
    return malloc(size);
}

static bool is_close(float a, float b, float rtol, float atol) {
  return fabs(a - b) <= (atol + rtol*fabs(b));
}

#endif
