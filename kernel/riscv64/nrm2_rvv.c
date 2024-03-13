/***************************************************************************
Copyright (c) 2022, The OpenBLAS Project
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

#if defined(DOUBLE)
#define VSETVL             __riscv_vsetvl_e64m4
#define FLOAT_V_T           vfloat64m4_t
#define FLOAT_V_T_M1        vfloat64m1_t
#define VLEV_FLOAT          __riscv_vle64_v_f64m4
#define VLSEV_FLOAT         __riscv_vlse64_v_f64m4
#define VFMVVF_FLOAT        __riscv_vfmv_v_f_f64m4
#define VFMVSF_FLOAT        __riscv_vfmv_s_f_f64m4
#define VFMVVF_FLOAT_M1     __riscv_vfmv_v_f_f64m1
#define MASK_T              vbool16_t
#define VFABS               __riscv_vfabs_v_f64m4
#define VMFNE               __riscv_vmfne_vf_f64m4_b16
#define VMFGT               __riscv_vmfgt_vv_f64m4_b16
#define VMFEQ               __riscv_vmfeq_vf_f64m4_b16
#define VCPOP               __riscv_vcpop_m_b16
#define VFREDMAX            __riscv_vfredmax_vs_f64m4_f64m1
#define VFREDMIN            __riscv_vfredmin_vs_f64m4_f64m1
#define VFIRST              __riscv_vfirst_m_b16
#define VRGATHER            __riscv_vrgather_vx_f64m4
#define VFDIV               __riscv_vfdiv_vv_f64m4
#define VFDIV_M             __riscv_vfdiv_vv_f64m4_mu
#define VFMUL               __riscv_vfmul_vv_f64m4
#define VFMUL_M             __riscv_vfmul_vv_f64m4_mu
#define VFMACC              __riscv_vfmacc_vv_f64m4
#define VFMACC_M            __riscv_vfmacc_vv_f64m4_mu
#define VMSBF               __riscv_vmsbf_m_b16
#define VMSOF               __riscv_vmsof_m_b16
#define VMAND               __riscv_vmand_mm_b16
#define VMANDN              __riscv_vmand_mm_b16
#define VFREDSUM            __riscv_vfredusum_vs_f64m4_f64m1
#define VMERGE              __riscv_vmerge_vvm_f64m4
#define VSEV_FLOAT          __riscv_vse64_v_f64m4
#define EXTRACT_FLOAT0_V(v) __riscv_vfmv_f_s_f64m4_f64(v)
#define ABS fabs
#else
#define VSETVL              __riscv_vsetvl_e32m4
#define FLOAT_V_T           vfloat32m4_t
#define FLOAT_V_T_M1        vfloat32m1_t
#define VLEV_FLOAT          __riscv_vle32_v_f32m4
#define VLSEV_FLOAT         __riscv_vlse32_v_f32m4
#define VFMVVF_FLOAT        __riscv_vfmv_v_f_f32m4
#define VFMVSF_FLOAT        __riscv_vfmv_s_f_f32m4
#define VFMVVF_FLOAT_M1     __riscv_vfmv_v_f_f32m1
#define MASK_T              vbool8_t
#define VFABS               __riscv_vfabs_v_f32m4
#define VMFNE               __riscv_vmfne_vf_f32m4_b8
#define VMFGT               __riscv_vmfgt_vv_f32m4_b8
#define VMFEQ               __riscv_vmfeq_vf_f32m4_b8
#define VCPOP               __riscv_vcpop_m_b8
#define VFREDMAX            __riscv_vfredmax_vs_f32m4_f32m1
#define VFREDMIN            __riscv_vfredmin_vs_f32m4_f32m1
#define VFIRST              __riscv_vfirst_m_b8
#define VRGATHER            __riscv_vrgather_vx_f32m4
#define VFDIV               __riscv_vfdiv_vv_f32m4
#define VFDIV_M             __riscv_vfdiv_vv_f32m4_mu
#define VFMUL               __riscv_vfmul_vv_f32m4
#define VFMUL_M             __riscv_vfmul_vv_f32m4_mu
#define VFMACC              __riscv_vfmacc_vv_f32m4
#define VFMACC_M            __riscv_vfmacc_vv_f32m4_mu
#define VMSBF               __riscv_vmsbf_m_b8
#define VMSOF               __riscv_vmsof_m_b8
#define VMAND               __riscv_vmand_mm_b8
#define VMANDN              __riscv_vmand_mm_b8
#define VFREDSUM            __riscv_vfredusum_vs_f32m4_f32m1
#define VMERGE              __riscv_vmerge_vvm_f32m4
#define VSEV_FLOAT          __riscv_vse32_v_f32m4
#define EXTRACT_FLOAT0_V(v) __riscv_vfmv_f_s_f32m4_f32(v)
#define ABS fabsf
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0;

	if (n <= 0 || inc_x == 0) return(0.0);
        if(n == 1) return (ABS(x[0]));

        unsigned int gvl = 0;

        MASK_T nonzero_mask;
        MASK_T scale_mask;

        gvl = VSETVL(n);
        FLOAT_V_T v0;
        FLOAT_V_T v_ssq = VFMVVF_FLOAT(0, gvl);
        FLOAT_V_T v_scale = VFMVVF_FLOAT(0, gvl);

        FLOAT scale = 0;
        FLOAT ssq = 0;
        unsigned int stride_x = inc_x * sizeof(FLOAT);
        int idx = 0;

        if( n >= gvl && inc_x > 0 ) // don't pay overheads if we're not doing useful work
        {
                for(i=0; i<n/gvl; i++){
                        v0 = VLSEV_FLOAT( &x[idx], stride_x, gvl );
                        nonzero_mask = VMFNE( v0, 0, gvl );
                        v0 = VFABS( v0, gvl );
                        scale_mask = VMFGT( v0, v_scale, gvl );

                        // assume scale changes are relatively infrequent

                        // unclear if the vcpop+branch is actually a win
                        // since the operations being skipped are predicated anyway
                        // need profiling to confirm
                        if( VCPOP(scale_mask, gvl) ) 
                        {
                                v_scale = VFDIV_M( scale_mask, v_scale, v_scale, v0, gvl );
                                v_scale = VFMUL_M( scale_mask, v_scale, v_scale, v_scale, gvl );
                                v_ssq = VFMUL_M( scale_mask, v_ssq, v_ssq, v_scale, gvl );
                                v_scale = VMERGE( v_scale, v0, scale_mask, gvl );
                        }
                        v0 = VFDIV_M( nonzero_mask, v0, v0, v_scale, gvl );
                        v_ssq = VFMACC_M( nonzero_mask, v_ssq, v0, v0, gvl );
                        idx += inc_x * gvl;
                }

                // we have gvl elements which we accumulated independently, with independent scales
                // we need to combine these
                // naive sort so we process small values first to avoid losing information
                // could use vector sort extensions where available, but we're dealing with gvl elts at most

                FLOAT * out_ssq = alloca(gvl*sizeof(FLOAT));
                FLOAT * out_scale = alloca(gvl*sizeof(FLOAT));
                VSEV_FLOAT( out_ssq, v_ssq, gvl );
                VSEV_FLOAT( out_scale, v_scale, gvl );
                for( int a = 0; a < (gvl-1); ++a )
                {
                        int smallest = a;
                        for( size_t b = a+1; b < gvl; ++b )
                                if( out_scale[b] < out_scale[smallest] )
                                        smallest = b;
                        if( smallest != a )
                        {
                                FLOAT tmp1 = out_ssq[a];
                                FLOAT tmp2 = out_scale[a];
                                out_ssq[a] = out_ssq[smallest];
                                out_scale[a] = out_scale[smallest];
                                out_ssq[smallest] = tmp1;
                                out_scale[smallest] = tmp2;
                        }
                }

                int a = 0;
                while( a<gvl && out_scale[a] == 0 )
                        ++a;

                if( a < gvl ) 
                {
                        ssq = out_ssq[a];
                        scale = out_scale[a];
                        ++a;
                        for( ; a < gvl; ++a ) 
                        {
                                ssq = ssq * ( scale / out_scale[a] ) * ( scale / out_scale[a] ) + out_ssq[a];
                                scale = out_scale[a];
                        }
                }
        }

        //finish any tail using scalar ops
        i*=gvl*inc_x;
        n*=inc_x;
        while(abs(i) < abs(n)){
                if ( x[i] != 0.0 ){
                        FLOAT absxi = ABS( x[i] );
                        if ( scale < absxi ){
                                ssq = 1 + ssq * ( scale / absxi ) * ( scale / absxi );
                                scale = absxi ;
                        }
                        else{
                                ssq += ( absxi/scale ) * ( absxi/scale );
                        }

                }

                i += inc_x;
        }

	return(scale * sqrt(ssq));
}


