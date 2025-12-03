/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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

#define M_BLOCKSIZE SGEMM_UNROLL_M
#define N_BLOCKSIZE SGEMM_UNROLL_N

#if M_BLOCKSIZE == 4
#define RVV_MUL __riscv_vfmul_vf_f32m1
#define RVV_MACC __riscv_vfmacc_vf_f32m1
#define RVV_LOAD __riscv_vle32_v_f32m1
#define RVV_STORE __riscv_vse32_v_f32m1
#define VECTOR_T vfloat32m1_t
#elif M_BLOCKSIZE == 8
#define RVV_MUL __riscv_vfmul_vf_f32m2
#define RVV_MACC __riscv_vfmacc_vf_f32m2
#define RVV_LOAD __riscv_vle32_v_f32m2
#define RVV_STORE __riscv_vse32_v_f32m2
#define VECTOR_T vfloat32m2_t
#elif M_BLOCKSIZE == 16
#define RVV_MUL __riscv_vfmul_vf_f32m4
#define RVV_MACC __riscv_vfmacc_vf_f32m4
#define RVV_LOAD __riscv_vle32_v_f32m4
#define RVV_STORE __riscv_vse32_v_f32m4
#define VECTOR_T vfloat32m4_t
#else
#error "Unsupported M_BLOCKSIZE value"
#endif

#include "gemm_kernel_rvv_vlv_common.h"
