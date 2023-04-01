/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
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

#if !defined(DOUBLE)
#define VSETVL(n)               __riscv_vsetvl_e32m4(n)
#define FLOAT_V_T               vfloat32m4_t
#define VLSEG_FLOAT             __riscv_vlseg2e32_v_f32m4
#define VSSEG_FLOAT             __riscv_vsseg2e32_v_f32m4
#define VFMVVF_FLOAT            __riscv_vfmv_v_f_f32m4
#define VFMULVF_FLOAT           __riscv_vfmul_vf_f32m4
#define VFADDVV_FLOAT           __riscv_vfadd_vv_f32m4
#define VFSUBVV_FLOAT           __riscv_vfsub_vv_f32m4
#else
#define VSETVL(n)               __riscv_vsetvl_e64m4(n)
#define FLOAT_V_T               vfloat64m4_t
#define VLSEG_FLOAT             __riscv_vlseg2e64_v_f64m4
#define VSSEG_FLOAT             __riscv_vsseg2e64_v_f64m4
#define VFMVVF_FLOAT            __riscv_vfmv_v_f_f64m4
#define VFMULVF_FLOAT           __riscv_vfmul_vf_f64m4
#define VFADDVV_FLOAT           __riscv_vfadd_vv_f64m4
#define VFSUBVV_FLOAT           __riscv_vfsub_vv_f64m4
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1,
          FLOAT beta_r, FLOAT beta_i,
          FLOAT *dummy2, BLASLONG dummy3,
          FLOAT *dummy4, BLASLONG dummy5,
          FLOAT *c, BLASLONG ldc)
{
    BLASLONG chunk;
    FLOAT *c_offset;
	size_t vl;
    FLOAT_V_T vr, vi, v1, v2, v3, v4;

    ldc *= 2;
    c_offset = c;

    if (beta_r == 0.0 && beta_i == 0.0) {

        vl = VSETVL(m);
        vr = VFMVVF_FLOAT(0.0, vl);
        vi = VFMVVF_FLOAT(0.0, vl);

        for( ; n > 0; n--, c += ldc) {
            c_offset = c;

            for(chunk=m; chunk > 0; chunk -= vl, c_offset += vl*2) {
                vl = VSETVL(chunk);

                VSSEG_FLOAT(c_offset, vr, vi, vl);
			}
		}

    } else {

        for( ; n > 0; n--, c += ldc) {
            c_offset = c;

            for(chunk=m; chunk > 0; chunk -= vl, c_offset += vl*2) {
                vl = VSETVL(chunk);

                VLSEG_FLOAT(&vr, &vi, c_offset, vl);

                v1 = VFMULVF_FLOAT(vr, beta_r, vl);
                v2 = VFMULVF_FLOAT(vi, beta_i, vl);

                v3 = VFMULVF_FLOAT(vi, beta_r, vl);
                v4 = VFMULVF_FLOAT(vr, beta_i, vl);

				vr = VFSUBVV_FLOAT(v1, v2, vl);
				vi = VFADDVV_FLOAT(v3, v4, vl);

                VSSEG_FLOAT(c_offset, vr, vi, vl);
			}
		}

	}

    return 0;
}
