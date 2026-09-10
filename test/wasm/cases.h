/*
Copyright (c) 2026, The OpenBLAS Project
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
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TEST_WASM_CASES_H
#define TEST_WASM_CASES_H

/* Dense remainders around 4x4 / 8x4 / 2x2 tiles, plus a few larger sizes. */
static const int SIZES_L1[] = {
    0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
    127, 128, 129, 255, 256, 257, 1023, 1024};
static const int NS_L1 = (int)(sizeof(SIZES_L1) / sizeof(SIZES_L1[0]));

static const int SIZES_L2[] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    31, 32, 33, 34, 35, 36, 63, 64, 65, 127, 128, 129};
static const int NS_L2 = (int)(sizeof(SIZES_L2) / sizeof(SIZES_L2[0]));

static const int SIZES_L3[] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    31, 32, 33, 34, 35, 36, 63, 64, 65, 127, 128, 129};
static const int NS_L3 = (int)(sizeof(SIZES_L3) / sizeof(SIZES_L3[0]));

/* Complex GEMM is heavier; keep max moderate while hitting 2x2 remainders. */
static const int SIZES_CZ[] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65};
static const int NS_CZ = (int)(sizeof(SIZES_CZ) / sizeof(SIZES_CZ[0]));

/* Full standard-BLAS coverage: broad remainder sampling without large cases. */
static const int SIZES_FULL[] = {
    1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32};
static const int NS_FULL = (int)(sizeof(SIZES_FULL) / sizeof(SIZES_FULL[0]));

static const int INCS[] = {0, 1, 2, 3};
static const int NINCS = (int)(sizeof(INCS) / sizeof(INCS[0]));

#endif /* TEST_WASM_CASES_H */
