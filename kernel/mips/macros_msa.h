/*******************************************************************************
Copyright (c) 2016, The OpenBLAS Project
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
*******************************************************************************/

#ifndef __MACROS_MSA_H__
#define __MACROS_MSA_H__

#include <msa.h>

#define LD_W(RTYPE, psrc) *((RTYPE *)(psrc))
#define LD_SP(...) LD_W(v4f32, __VA_ARGS__)

#define LD_D(RTYPE, psrc) *((RTYPE *)(psrc))
#define LD_DP(...) LD_D(v2f64, __VA_ARGS__)

#define ST_W(RTYPE, in, pdst) *((RTYPE *)(pdst)) = (in)
#define ST_SP(...) ST_W(v4f32, __VA_ARGS__)

#define ST_D(RTYPE, in, pdst) *((RTYPE *)(pdst)) = (in)
#define ST_DP(...) ST_D(v2f64, __VA_ARGS__)

/* Description : Load 2 vectors of single precision floating point elements with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1
                 Return Type - single precision floating point
*/
#define LD_SP2(psrc, stride, out0, out1)  \
{                                         \
    out0 = LD_SP((psrc));                 \
    out1 = LD_SP((psrc) + stride);        \
}

/* Description : Load 2 vectors of double precision floating point elements with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1
                 Return Type - double precision floating point
*/
#define LD_DP2(psrc, stride, out0, out1)  \
{                                         \
    out0 = LD_DP((psrc));                 \
    out1 = LD_DP((psrc) + stride);        \
}

#define LD_DP4(psrc, stride, out0, out1, out2, out3)  \
{                                                     \
    LD_DP2(psrc, stride, out0, out1)                  \
    LD_DP2(psrc + 2 * stride, stride, out2, out3)     \
}

/* Description : Store vectors of single precision floating point elements with stride
   Arguments   : Inputs - in0, in1, pdst, stride
   Details     : Store 4 single precision floating point elements from 'in0' to (pdst)
                 Store 4 single precision floating point elements from 'in1' to (pdst + stride)
*/
#define ST_SP2(in0, in1, pdst, stride)  \
{                                       \
    ST_SP(in0, (pdst));                 \
    ST_SP(in1, (pdst) + stride);        \
}

#define ST_SP4(in0, in1, in2, in3, pdst, stride)    \
{                                                   \
    ST_SP2(in0, in1, (pdst), stride);               \
    ST_SP2(in2, in3, (pdst + 2 * stride), stride);  \
}

#define ST_SP8(in0, in1, in2, in3, in4, in5, in6, in7, pdst, stride)  \
{                                                                     \
    ST_SP4(in0, in1, in2, in3, (pdst), stride);                       \
    ST_SP4(in4, in5, in6, in7, (pdst + 4 * stride), stride);          \
}

/* Description : Store vectors of double precision floating point elements with stride
   Arguments   : Inputs - in0, in1, pdst, stride
   Details     : Store 2 double precision floating point elements from 'in0' to (pdst)
                 Store 2 double precision floating point elements from 'in1' to (pdst + stride)
*/
#define ST_DP2(in0, in1, pdst, stride)  \
{                                       \
    ST_DP(in0, (pdst));                 \
    ST_DP(in1, (pdst) + stride);        \
}

#define ST_DP4(in0, in1, in2, in3, pdst, stride)   \
{                                                  \
    ST_DP2(in0, in1, (pdst), stride);              \
    ST_DP2(in2, in3, (pdst) + 2 * stride, stride); \
}

#define ST_DP8(in0, in1, in2, in3, in4, in5, in6, in7, pdst, stride)  \
{                                                                     \
    ST_DP4(in0, in1, in2, in3, (pdst), stride);                       \
    ST_DP4(in4, in5, in6, in7, (pdst) + 4 * stride, stride);          \
}

/* Description : Interleave both left and right half of input vectors
   Arguments   : Inputs  - in0, in1
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Right half of byte elements from 'in0' and 'in1' are
                 interleaved and written to 'out0'
*/
#define ILVRL_W2(RTYPE, in0, in1, out0, out1)               \
{                                                           \
    out0 = (RTYPE) __msa_ilvr_w((v4i32) in0, (v4i32) in1);  \
    out1 = (RTYPE) __msa_ilvl_w((v4i32) in0, (v4i32) in1);  \
}
#define ILVRL_W2_SW(...) ILVRL_W2(v4i32, __VA_ARGS__)

#define ILVRL_D2(RTYPE, in0, in1, out0, out1)               \
{                                                           \
    out0 = (RTYPE) __msa_ilvr_d((v2i64) in0, (v2i64) in1);  \
    out1 = (RTYPE) __msa_ilvl_d((v2i64) in0, (v2i64) in1);  \
}
#define ILVRL_D2_DP(...) ILVRL_D2(v2f64, __VA_ARGS__)

/* Description : Indexed word element values are replicated to all
                 elements in output vector
   Arguments   : Inputs  - in, stidx
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : 'stidx' element value from 'in' vector is replicated to all
                 elements in 'out0' vector
                 'stidx + 1' element value from 'in' vector is replicated to all
                 elements in 'out1' vector
                 Valid index range for word operation is 0-3
*/
#define SPLATI_W2(RTYPE, in, stidx, out0, out1)            \
{                                                          \
    out0 = (RTYPE) __msa_splati_w((v4i32) in, stidx);      \
    out1 = (RTYPE) __msa_splati_w((v4i32) in, (stidx+1));  \
}

#define SPLATI_W4(RTYPE, in, out0, out1, out2, out3)  \
{                                                     \
    SPLATI_W2(RTYPE, in, 0, out0, out1);              \
    SPLATI_W2(RTYPE, in, 2, out2, out3);              \
}
#define SPLATI_W4_SP(...) SPLATI_W4(v4f32, __VA_ARGS__)

/* Description : Transpose 4x4 block with word elements in vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1, out2, out3
                 Return Type - as per RTYPE
*/
#define TRANSPOSE4x4_W(RTYPE, in0, in1, in2, in3, out0, out1, out2, out3)  \
{                                                                   \
    v4i32 s0_m, s1_m, s2_m, s3_m;                                   \
                                                                    \
    ILVRL_W2_SW(in1, in0, s0_m, s1_m);                              \
    ILVRL_W2_SW(in3, in2, s2_m, s3_m);                              \
                                                                    \
    out0 = (RTYPE) __msa_ilvr_d((v2i64) s2_m, (v2i64) s0_m);        \
    out1 = (RTYPE) __msa_ilvl_d((v2i64) s2_m, (v2i64) s0_m);        \
    out2 = (RTYPE) __msa_ilvr_d((v2i64) s3_m, (v2i64) s1_m);        \
    out3 = (RTYPE) __msa_ilvl_d((v2i64) s3_m, (v2i64) s1_m);        \
}

#define TRANSPOSE4x4_SP_SP(...) TRANSPOSE4x4_W(v4f32, __VA_ARGS__)

#endif  /* __MACROS_MSA_H__ */
