/**********************************Zero Vectors**************************************************/

.macro ZERO_CVEC_8x4
    vzero  %v16  
    vzero  %v17   
    vzero  %v18   
    vzero  %v19   
    vzero  %v20   
    vzero  %v21   
    vzero  %v22   
    vzero  %v23   
    vzero  %v24   
    vzero  %v25   
    vzero  %v26   
    vzero  %v27   
    vzero  %v28   
    vzero  %v29   
    vzero  %v30   
    vzero  %v31 
.endm


.macro ZERO_CVEC_8x2
    vzero  %v16  
    vzero  %v17   
    vzero  %v18   
    vzero  %v19   
    vzero  %v20   
    vzero  %v21   
    vzero  %v22   
    vzero  %v23   
     
.endm

.macro ZERO_CVEC_8x1
    vzero  %v16  
    vzero  %v17
    vzero  %v18  
    vzero  %v19     
.endm

.macro ZERO_CVEC_4x4
    vzero  %v16  
    vzero  %v17    
    vzero  %v20   
    vzero  %v21    
    vzero  %v24   
    vzero  %v25     
    vzero  %v28   
    vzero  %v29   
.endm

.macro ZERO_CVEC_4x2
    vzero  %v16  
    vzero  %v17     
    vzero  %v20   
    vzero  %v21    
     
.endm

.macro ZERO_CVEC_4x1
    vzero  %v16  
    vzero  %v17      
.endm

.macro ZERO_CVEC_2x4
    vzero  %v16  
    vzero  %v17     
    vzero  %v20   
    vzero  %v21    
     
.endm

.macro ZERO_CVEC_2x2
    vzero  %v16       
    vzero  %v20      
     
.endm

.macro ZERO_CVEC_2x1
    vzero  %v16   
.endm

.macro ZERO_CVEC_1x4
    vzero  %v16  
    vzero  %v17      
.endm

.macro ZERO_CVEC_1x2
    vzero  %v16       
.endm

.macro ZERO_CVEC_1x1
    LZDR %f1       
.endm


/***********************************Helper Calculations*************************************/
#define unit_size 8
#define DISP(ind,stride,disp) (ind*stride+disp)
#define DISP8(ind,disp) (ind*unit_size*8+disp)
#define DISP4(ind,disp) (ind*unit_size*4+disp)
#define DISP2(ind,disp) (ind*unit_size*2+disp)
#define DISP1(ind,disp) (ind*unit_size+disp)
#define N8  (8*unit_size)
#define N4  (4*unit_size)
#define N2  (2*unit_size)
#define N1  (1*unit_size)

.macro Calculate_8x4_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP4(\Index ,0)(\PTR_B_REG)  
        vlrepg %v1, DISP4(\Index ,8)(\PTR_B_REG)    
        vl  %v2,    DISP8(\Index , 0)(\PTR_A_REG)     
        vl  %v3,    DISP8(\Index ,16)(\PTR_A_REG)     
        vl  %v4,    DISP8(\Index ,32)(\PTR_A_REG)   
        vl  %v5,    DISP8(\Index ,48)(\PTR_A_REG)      
        vfmadb    %v16,%v2,%v7,%v16   
        vfmadb    %v17,%v3,%v7,%v17   
        vfmadb    %v18,%v4,%v7,%v18  
        vfmadb    %v19,%v5,%v7,%v19     
        vfmadb    %v20,%v2,%v1,%v20   
        vfmadb    %v21,%v3,%v1,%v21 
        vfmadb    %v22,%v4,%v1,%v22  
        vfmadb    %v23,%v5,%v1,%v23  
        vlrepg %v7, DISP4(\Index ,16)(\PTR_B_REG)      
        vlrepg %v1, DISP4(\Index ,24)(\PTR_B_REG) 
   .if \IsLast==1         
        la \PTR_A_REG, DISP8(\Index ,64)(\PTR_A_REG)   
   .endif    
        vfmadb    %v24,%v2,%v7,%v24   
        vfmadb    %v25,%v3,%v7,%v25   
        vfmadb    %v26,%v4,%v7,%v26  
        vfmadb    %v27,%v5,%v7,%v27        
        vfmadb    %v28,%v2,%v1,%v28   
        vfmadb    %v29,%v3,%v1,%v29   
        vfmadb    %v30,%v4,%v1,%v30  
        vfmadb    %v31,%v5,%v1,%v31
    .if \IsLast==1      
        la \PTR_B_REG, DISP4(\Index ,32)(\PTR_B_REG)
    .endif    
.endm

.macro Calculate_8x2_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP2(\Index ,0)(\PTR_B_REG)  
        vlrepg %v1, DISP2(\Index ,8)(\PTR_B_REG)    
        vl  %v2,    DISP8(\Index ,0)(\PTR_A_REG)     
        vl  %v3,    DISP8(\Index ,16)(\PTR_A_REG)     
        vl  %v4,    DISP8(\Index ,32)(\PTR_A_REG)   
        vl  %v5,    DISP8(\Index ,48)(\PTR_A_REG)      
        vfmadb    %v16,%v2,%v7,%v16   
        vfmadb    %v17,%v3,%v7,%v17   
        vfmadb    %v18,%v4,%v7,%v18  
        vfmadb    %v19,%v5,%v7,%v19     
        vfmadb    %v20,%v2,%v1,%v20  
        vfmadb    %v21,%v3,%v1,%v21
    .if \IsLast==1             
        la \PTR_A_REG, DISP8(\Index ,64)(\PTR_A_REG)  
    .endif        
        vfmadb    %v22,%v4,%v1,%v22  
        vfmadb    %v23,%v5,%v1,%v23 
    .if \IsLast==1                
        la \PTR_B_REG, DISP2(\Index ,16)(\PTR_B_REG)
    .endif         
.endm

.macro Calculate_8x1_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP1(\Index ,0)(\PTR_B_REG)    
        vl  %v2,    DISP8(\Index ,0)(\PTR_A_REG)     
        vl  %v3,    DISP8(\Index ,16)(\PTR_A_REG)     
        vl  %v4,    DISP8(\Index ,32)(\PTR_A_REG)   
        vl  %v5,    DISP8(\Index ,48)(\PTR_A_REG)     
        vfmadb    %v16,%v2,%v7,%v16      
    .if \IsLast==1            
        la \PTR_B_REG, DISP1(\Index ,8)(\PTR_B_REG) 
    .endif    
        vfmadb    %v17,%v3,%v7,%v17   
        vfmadb    %v18,%v4,%v7,%v18  
        vfmadb    %v19,%v5,%v7,%v19    
    .if \IsLast==1    
        la \PTR_A_REG, DISP8(\Index ,64)(\PTR_A_REG)
    .endif      
.endm

.macro Calculate_4x4_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP4(\Index ,0)(\PTR_B_REG)  
        vlrepg %v1, DISP4(\Index ,8)(\PTR_B_REG)    
        vl  %v2,    DISP4(\Index ,0)(\PTR_A_REG)     
        vl  %v3,    DISP4(\Index ,16)(\PTR_A_REG)         
        vfmadb    %v16,%v2,%v7,%v16   
        vfmadb    %v17,%v3,%v7,%v17       
        vfmadb    %v20,%v2,%v1,%v20   
        vfmadb    %v21,%v3,%v1,%v21 
        vlrepg %v7, DISP4(\Index ,16)(\PTR_B_REG)      
        vlrepg %v1, DISP4(\Index ,24)(\PTR_B_REG)  
    .if \IsLast==1           
        la \PTR_A_REG, DISP4(\Index ,32)(\PTR_A_REG)
    .endif           
        vfmadb    %v24,%v2,%v7,%v24   
        vfmadb    %v25,%v3,%v7,%v25  
        vfmadb    %v28,%v2,%v1,%v28   
        vfmadb    %v29,%v3,%v1,%v29 
    .if \IsLast==1               
        la \PTR_B_REG, DISP4(\Index ,32)(\PTR_B_REG) 
    .endif      
.endm

.macro Calculate_4x2_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP2(\Index ,0)(\PTR_B_REG)  
        vlrepg %v1, DISP2(\Index ,8)(\PTR_B_REG)    
        vl  %v2, DISP4(\Index ,0)(\PTR_A_REG)     
        vl  %v3, DISP4(\Index ,16)(\PTR_A_REG)       
        vfmadb    %v16,%v2,%v7,%v16   
        vfmadb    %v17,%v3,%v7,%v17 
    .if \IsLast==1    
        la \PTR_B_REG, DISP2(\Index ,16)(\PTR_B_REG)  
    .endif          
        vfmadb    %v20,%v2,%v1,%v20   
        vfmadb    %v21,%v3,%v1,%v21  
    .if \IsLast==1       
        la \PTR_A_REG, DISP4(\Index ,32)(\PTR_A_REG)   
    .endif      
.endm

.macro Calculate_4x1_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP1(\Index ,0)(\PTR_B_REG)    
        vl  %v2, DISP4(\Index ,0)(\PTR_A_REG)     
        vl  %v3, DISP4(\Index ,16)(\PTR_A_REG) 
    .if \IsLast==1                 
        la \PTR_B_REG, DISP1(\Index ,8)(\PTR_B_REG)  
    .endif      
        vfmadb    %v16,%v2,%v7,%v16      
        vfmadb    %v17,%v3,%v7,%v17   
    .if \IsLast==1          
        la \PTR_A_REG, DISP4(\Index ,32)(\PTR_A_REG)
    .endif         
.endm

.macro Calculate_2x2_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP2(\Index ,0)(\PTR_B_REG)  
        vlrepg %v1, DISP2(\Index ,8)(\PTR_B_REG)    
        vl  %v2,    DISP2(\Index ,0)(\PTR_A_REG)         
        vfmadb    %v16,%v2,%v7,%v16    
    .if \IsLast==1           
        la \PTR_A_REG, DISP2(\Index ,16)(\PTR_A_REG)  
    .endif         
        vfmadb    %v20,%v2,%v1,%v20  
    .if \IsLast==1          
        la \PTR_B_REG, DISP2(\Index ,16)(\PTR_B_REG)
    .endif        
.endm
 
 

.macro Calculate_2x1_I  PTR_A_REG,PTR_B_REG,Index,IsLast
        vlrepg %v7, DISP1(\Index ,0)(\PTR_B_REG)    
        vl  %v2, DISP2(\Index ,0)(\PTR_A_REG)  
    .if \IsLast==1              
        la \PTR_B_REG, DISP1(\Index ,8)(\PTR_B_REG)
    .endif        
        vfmadb    %v16,%v2,%v7,%v16 
    .if \IsLast==1     
        la \PTR_A_REG, DISP2(\Index ,16)(\PTR_A_REG)
    .endif         
.endm

.macro  Calculate_1x1_I  PTR_A_REG,PTR_B_REG,Index,IsLast
     ld %f2,DISP1(\Index ,0)(\PTR_A_REG) /**a*/
    .if \IsLast==1  
      la  \PTR_A_REG,DISP1(\Index ,8)(\PTR_A_REG)
    .endif 
     madb %f1,%f2,DISP1(\Index ,0)(\PTR_B_REG)
    .if \IsLast==1       
      la  \PTR_B_REG,DISP1(\Index ,8)(\PTR_B_REG) 
    .endif    
.endm

.macro CALC_8x4  PTR_A_REG,PTR_B_REG
   Calculate_8x4_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_8x4_4  PTR_A_REG,PTR_B_REG 
   Calculate_8x4_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_8x4_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_8x4_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_8x4_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_8x2  PTR_A_REG,PTR_B_REG
   Calculate_8x2_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_8x2_4  PTR_A_REG,PTR_B_REG 
   Calculate_8x2_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_8x2_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_8x2_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_8x2_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_8x1  PTR_A_REG,PTR_B_REG
   Calculate_8x1_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_8x1_4  PTR_A_REG,PTR_B_REG 
   Calculate_8x1_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_8x1_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_8x1_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_8x1_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_4x4  PTR_A_REG,PTR_B_REG
   Calculate_4x4_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_4x4_4  PTR_A_REG,PTR_B_REG 
   Calculate_4x4_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_4x4_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_4x4_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_4x4_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_4x2  PTR_A_REG,PTR_B_REG
   Calculate_4x2_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_4x2_4  PTR_A_REG,PTR_B_REG 
   Calculate_4x2_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_4x2_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_4x2_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_4x2_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_4x1  PTR_A_REG,PTR_B_REG
   Calculate_4x1_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_4x1_4  PTR_A_REG,PTR_B_REG 
   Calculate_4x1_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_4x1_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_4x1_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_4x1_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_2x4  PTR_A_REG,PTR_B_REG
   Calculate_4x2_I  \PTR_B_REG,\PTR_A_REG,0,1
.endm

.macro CALC_2x4_4  PTR_A_REG,PTR_B_REG 
   Calculate_4x2_I  \PTR_B_REG,\PTR_A_REG,0,0
   Calculate_4x2_I  \PTR_B_REG,\PTR_A_REG,1,0
   Calculate_4x2_I  \PTR_B_REG,\PTR_A_REG,2,0
   Calculate_4x2_I  \PTR_B_REG,\PTR_A_REG,3,1 
.endm

.macro CALC_2x2  PTR_A_REG,PTR_B_REG
   Calculate_2x2_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_2x2_4  PTR_A_REG,PTR_B_REG 
   Calculate_2x2_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_2x2_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_2x2_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_2x2_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_2x1  PTR_A_REG,PTR_B_REG
   Calculate_2x1_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_2x1_4  PTR_A_REG,PTR_B_REG 
   Calculate_2x1_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_2x1_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_2x1_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_2x1_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm

.macro CALC_1x4  PTR_A_REG,PTR_B_REG
   Calculate_4x1_I  \PTR_B_REG,\PTR_A_REG,0,1
.endm

.macro CALC_1x4_4  PTR_A_REG,PTR_B_REG 
   Calculate_4x1_I  \PTR_B_REG,\PTR_A_REG,0,0
   Calculate_4x1_I  \PTR_B_REG,\PTR_A_REG,1,0
   Calculate_4x1_I  \PTR_B_REG,\PTR_A_REG,2,0
   Calculate_4x1_I  \PTR_B_REG,\PTR_A_REG,3,1 
.endm

.macro CALC_1x2  PTR_A_REG,PTR_B_REG
   Calculate_2x1_I  \PTR_B_REG,\PTR_A_REG,0,1
.endm

.macro CALC_1x2_4  PTR_A_REG,PTR_B_REG 
   Calculate_2x1_I  \PTR_B_REG,\PTR_A_REG,0,0
   Calculate_2x1_I  \PTR_B_REG,\PTR_A_REG,1,0
   Calculate_2x1_I  \PTR_B_REG,\PTR_A_REG,2,0
   Calculate_2x1_I  \PTR_B_REG,\PTR_A_REG,3,1 
.endm

.macro CALC_1x1  PTR_A_REG,PTR_B_REG
   Calculate_1x1_I  \PTR_A_REG,\PTR_B_REG,0,1
.endm

.macro CALC_1x1_4  PTR_A_REG,PTR_B_REG 
   Calculate_1x1_I  \PTR_A_REG,\PTR_B_REG,0,0
   Calculate_1x1_I  \PTR_A_REG,\PTR_B_REG,1,0
   Calculate_1x1_I  \PTR_A_REG,\PTR_B_REG,2,0
   Calculate_1x1_I  \PTR_A_REG,\PTR_B_REG,3,1 
.endm
 

/**************************************STORAGE*************************************************/


.macro Multiply_8x1 vr1,vr2,vr3,vr4,va1,va2,va3,va4,vb1   
 #if defined(TRMMKERNEL)
        vfmdb    \vr1,\va1,\vb1    
        vfmdb    \vr2,\va2,\vb1    
        vfmdb    \vr3,\va3,\vb1  
        vfmdb    \vr4,\va4,\vb1   
 #else
        vfmadb    \vr1,\va1,\vb1,\vr1    
        vfmadb    \vr2,\va2,\vb1,\vr2   
        vfmadb    \vr3,\va3,\vb1,\vr3  
        vfmadb    \vr4,\va4,\vb1,\vr4   
 #endif         
.endm    

.macro Multiply_4x1 vr1,vr2, va1,va2, vb1   
 #if defined(TRMMKERNEL)
        vfmdb    \vr1,\va1,\vb1    
        vfmdb    \vr2,\va2,\vb1     
 #else
        vfmadb    \vr1,\va1,\vb1,\vr1    
        vfmadb    \vr2,\va2,\vb1,\vr2    
 #endif        
.endm

.macro Multiply_2x1 vr1, va1,vb1   
 #if defined(TRMMKERNEL)
        vfmdb    \vr1,\va1,\vb1     
 #else
        vfmadb    \vr1,\va1,\vb1,\vr1      
 #endif         
.endm   

 
.macro STORE_8x4  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL , LV1 ,LV2
    la \LV1,0(\LDC_BYTE_ORIGINAL, \LDC_BYTE_ORIGINAL)
#if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG)
    vl %v2,16(\CIJ_REG)
    vl %v3,32(\CIJ_REG)
    vl %v4,48(\CIJ_REG)
#endif    
    Multiply_8x1 %v1,%v2,%v3,%v4,  %v16,%v17,%v18,%v19  ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG) 
    vst %v2,16(\CIJ_REG) 
    vst %v3,32(\CIJ_REG) 
    vst %v4,48(\CIJ_REG) 
    
    la \LV2,0(\LV1,\LDC_BYTE_ORIGINAL )
#if !defined(TRMMKERNEL)    
    vl %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v18,32(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v19,48(\CIJ_REG,\LDC_BYTE_ORIGINAL)
#endif
    Multiply_8x1 %v16,%v17,%v18,%v19,  %v20,%v21,%v22,%v23  ,\ALPHA_VECREG
    vst %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v18,32(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v19,48(\CIJ_REG,\LDC_BYTE_ORIGINAL)

 #if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG,\LV1)
    vl %v2,16(\CIJ_REG,\LV1)
    vl %v3,32(\CIJ_REG,\LV1)
    vl %v4,48(\CIJ_REG,\LV1)
#endif    
    Multiply_8x1 %v1,%v2,%v3,%v4,  %v24,%v25,%v26,%v27  ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG,\LV1) 
    vst %v2,16(\CIJ_REG,\LV1) 
    vst %v3,32(\CIJ_REG,\LV1) 
    vst %v4,48(\CIJ_REG,\LV1)  

#if !defined(TRMMKERNEL)    
    vl %v16,0(\CIJ_REG,\LV2)
    vl %v17,16(\CIJ_REG,\LV2)
    vl %v18,32(\CIJ_REG,\LV2)
    vl %v19,48(\CIJ_REG,\LV2)
#endif
    Multiply_8x1 %v16,%v17,%v18,%v19,  %v28,%v29,%v30,%v31  ,\ALPHA_VECREG
    vst %v16,0(\CIJ_REG,\LV2) 
    vst %v17,16(\CIJ_REG,\LV2) 
    vst %v18,32(\CIJ_REG,\LV2) 
    vst %v19,48(\CIJ_REG,\LV2) 

    la \CIJ_REG,64(\CIJ_REG)
 
.endm

.macro STORE_8x2  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL
#if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG)
    vl %v2,16(\CIJ_REG)
    vl %v3,32(\CIJ_REG)
    vl %v4,48(\CIJ_REG)
#endif    
    Multiply_8x1 %v1,%v2,%v3,%v4,  %v16,%v17,%v18,%v19  ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG) 
    vst %v2,16(\CIJ_REG) 
    vst %v3,32(\CIJ_REG) 
    vst %v4,48(\CIJ_REG) 
     
#if !defined(TRMMKERNEL)    
    vl %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v18,32(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v19,48(\CIJ_REG,\LDC_BYTE_ORIGINAL)
#endif
    Multiply_8x1 %v16,%v17,%v18,%v19,  %v20,%v21,%v22,%v23  ,\ALPHA_VECREG
    vst %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v18,32(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v19,48(\CIJ_REG,\LDC_BYTE_ORIGINAL)
 
   la \CIJ_REG,64(\CIJ_REG)
 
.endm

.macro STORE_8x1  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL
    
#if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG)
    vl %v2,16(\CIJ_REG)
    vl %v3,32(\CIJ_REG)
    vl %v4,48(\CIJ_REG)
#endif    
    Multiply_8x1 %v1,%v2,%v3,%v4,  %v16,%v17,%v18,%v19  ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG) 
    vst %v2,16(\CIJ_REG) 
    vst %v3,32(\CIJ_REG) 
    vst %v4,48(\CIJ_REG) 

    la \CIJ_REG,64(\CIJ_REG)
.endm


.macro STORE_4x4  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL,   LV1 ,LV2
    la \LV1,0(\LDC_BYTE_ORIGINAL, \LDC_BYTE_ORIGINAL)
#if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG)
    vl %v2,16(\CIJ_REG) 
#endif    
    Multiply_4x1 %v1,%v2 ,  %v16,%v17   ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG) 
    vst %v2,16(\CIJ_REG)  
    
    la \LV2,0(\LV1,\LDC_BYTE_ORIGINAL )
#if !defined(TRMMKERNEL)    
    vl %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
#endif
    Multiply_4x1 %v16,%v17 ,  %v20,%v21  ,\ALPHA_VECREG
    vst %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL)  

 #if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG,\LV1)
    vl %v2,16(\CIJ_REG,\LV1) 
#endif    
    Multiply_4x1 %v1,%v2 ,  %v24,%v25 ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG,\LV1) 
    vst %v2,16(\CIJ_REG,\LV1)    

#if !defined(TRMMKERNEL)    
    vl %v16,0(\CIJ_REG,\LV2)
    vl %v17,16(\CIJ_REG,\LV2) 
#endif
    Multiply_4x1 %v16,%v17,   %v28,%v29   ,\ALPHA_VECREG
    vst %v16,0(\CIJ_REG,\LV2) 
    vst %v17,16(\CIJ_REG,\LV2)   

   la \CIJ_REG,32(\CIJ_REG)
 
.endm


.macro STORE_4x2  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL
    
#if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG)
    vl %v2,16(\CIJ_REG) 
#endif    
    Multiply_4x1 %v1,%v2 ,  %v16,%v17   ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG) 
    vst %v2,16(\CIJ_REG)   

#if !defined(TRMMKERNEL)    
    vl %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    vl %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
#endif
    Multiply_4x1 %v16,%v17 ,  %v20,%v21  ,\ALPHA_VECREG
    vst %v16,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
    vst %v17,16(\CIJ_REG,\LDC_BYTE_ORIGINAL)   

    la \CIJ_REG,32(\CIJ_REG)
 
.endm
.macro STORE_4x1  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL
    
#if !defined(TRMMKERNEL)
    vl %v1,0(\CIJ_REG)
    vl %v2,16(\CIJ_REG) 
#endif    
    Multiply_4x1 %v1,%v2 ,  %v16,%v17   ,\ALPHA_VECREG
    vst %v1,0(\CIJ_REG) 
    vst %v2,16(\CIJ_REG)   
 
      la \CIJ_REG,32(\CIJ_REG)
 
.endm

.macro STORE_2x2  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL
 
#if !defined(TRMMKERNEL)    
    vl %v1,0(\CIJ_REG)
#endif    
    Multiply_2x1   %v1,%v16,\ALPHA_VECREG 
    vst %v1,0(\CIJ_REG) 

#if !defined(TRMMKERNEL) 
    vl %v2,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
#endif    
    Multiply_2x1   %v2,%v20,\ALPHA_VECREG 
    vst %v2,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 

    la \CIJ_REG,16(\CIJ_REG)
 
.endm


.macro STORE_2x1  ALPHA_VECREG,CIJ_REG , LDC_BYTE_ORIGINAL
 
#if !defined(TRMMKERNEL)    
    vl %v1,0(\CIJ_REG)
#endif    
    Multiply_2x1   %v1,%v16,\ALPHA_VECREG 
    vst %v1,0(\CIJ_REG) 

    la \CIJ_REG,16(\CIJ_REG)
.endm


/*STORE C1X1*/
.macro STORE_1x1  ALPHA_FLOAT,CIJ_REG,LDC_BYTE_ORIGINAL  

#if  defined(TRMMKERNEL) 
     mdbr %f1,\ALPHA_FLOAT
     std %f1,0(CIJ_LOCAL)  
#else    
     ld %f2,0(CIJ_LOCAL)
     madbr %f2,%f1,\ALPHA_FLOAT
     std %f2,0(CIJ_LOCAL) 
#endif     
     la \CIJ_REG,8(\CIJ_REG) 
.endm

/*reversed ones*/

.macro STORE_2x4  ALPHA_REG,CIJ_REG , LDC_BYTE_ORIGINAL , LV1 ,LV2
/**/
    vfmdb %v1,%v16,\ALPHA_REG
    vfmdb %v2,%v17,\ALPHA_REG
    vfmdb %v6,%v20,\ALPHA_REG
    vfmdb %v7,%v21,\ALPHA_REG
    vrepg %v4,%v1,1
    vrepg %v5,%v6,1
    la \LV1,0(\LDC_BYTE_ORIGINAL, \LDC_BYTE_ORIGINAL) 
 #if !defined(TRMMKERNEL)     
    adb %f1, 0(\CIJ_REG)
 #endif   
    std %f1,0(\CIJ_REG)
 #if !defined(TRMMKERNEL) 
    adb %f6, 8(\CIJ_REG)
 #endif   
    std %f6,8(\CIJ_REG)
  #if !defined(TRMMKERNEL)     
    adb %f4,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
  #endif  
    std %f4,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
  #if !defined(TRMMKERNEL)   
    adb %f5,8(\CIJ_REG,\LDC_BYTE_ORIGINAL)
  #endif  
    std %f5,8(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
     
      /*add LDC_BYTE */
    la \LV2,0(\LV1,\LDC_BYTE_ORIGINAL )
    vrepg %v4,%v2,1
    vrepg %v5,%v7,1
  #if !defined(TRMMKERNEL) 
    adb %f2,0(\CIJ_REG,\LV1)
  #endif  
    std %f2,0(\CIJ_REG,\LV1)
  #if !defined(TRMMKERNEL)   
    adb %f7,8(\CIJ_REG,\LV1)
  #endif  
    std %f7,8(\CIJ_REG,\LV1)
  #if  !defined(TRMMKERNEL) 
    adb %f4,0(\CIJ_REG,\LV2)
  #endif  
    std %f4,0(\CIJ_REG,\LV2)
  #if !defined(TRMMKERNEL) 
    adb %f5,8(\CIJ_REG,\LV2)
  #endif  
    std %f5,8(\CIJ_REG,\LV2)
    
     la \CIJ_REG,16(\CIJ_REG)

.endm

.macro STORE_1x4  ALPHA_REG,CIJ_REG , LDC_BYTE_ORIGINAL  , LV1 ,LV2
 
    vfmdb %v1,%v16,\ALPHA_REG
    vfmdb %v2,%v17,\ALPHA_REG
    vrepg %v4,%v1,1
    vrepg %v5,%v2,1
    la \LV1,0(\LDC_BYTE_ORIGINAL, \LDC_BYTE_ORIGINAL) 
  #if !defined(TRMMKERNEL)    
    adb %f1, 0(\CIJ_REG)
  #endif   
    std %f1,0(\CIJ_REG)
  #if !defined(TRMMKERNEL)     
    adb %f4,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
  #endif  
    std %f4,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
      /*add LDC_BYTE */
    la \LV2,0(\LV1,\LDC_BYTE_ORIGINAL )
  #if !defined(TRMMKERNEL)  
    adb %f2,0(\CIJ_REG,\LV1)
  #endif  
    std %f2,0(\CIJ_REG,\LV1)
  #if !defined(TRMMKERNEL)  
    adb %f5,0(\CIJ_REG,\LV2)
  #endif  
    std %f5,0(\CIJ_REG,\LV2)

    la \CIJ_REG,8(\CIJ_REG)

.endm

 .macro STORE_1x2  ALPHA_REG,CIJ_REG , LDC_BYTE_ORIGINAL
/**/
    vfmdb %v1,%v16,\ALPHA_REG 
    vrepg %v4,%v1,1  
   #if !defined(TRMMKERNEL)     
    adb %f1, 0(\CIJ_REG)
   #endif  
    std %f1,0(\CIJ_REG)
    #if !defined(TRMMKERNEL)  
    adb %f4,0(\CIJ_REG,\LDC_BYTE_ORIGINAL)
    #endif
    std %f4,0(\CIJ_REG,\LDC_BYTE_ORIGINAL) 
      
    la \CIJ_REG,8(\CIJ_REG)

.endm




/****************************TRMM POINTER REFRESH MACROSES*************************/

.macro RefreshPointers  PTR_A,PTR_B,OFF_VAL,B_VAL,C_A,C_B
    #if (defined(LEFT) &&  defined(TRANSA)) ||  (!defined(LEFT) && !defined(TRANSA))
               /*  ptrbb = bb;*/
            lgr \PTR_B,\B_VAL      /*refresh BPOINT*/

    #else
             /*  ptrba  =ptrba+ off*C_A;
              ptrbb = bb + off*C_B;*/
.if \C_B==4 
		.if \C_A==8 
			sllg \PTR_B, \OFF_VAL,5
			la \PTR_A,0(\PTR_A,\PTR_B) /*ptrba+off*4*/
			agr \PTR_A,\PTR_B /*ptrba+off*4**/
			la \PTR_B,0(\B_VAL,\PTR_B)	  
		.elseif \C_A==4
			sllg \PTR_B, \OFF_VAL,5 
			agr \PTR_A,\PTR_B /*ptrba+off*4**/
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.elseif \C_A==2
			sllg \PTR_B, \OFF_VAL,4 
			la \PTR_A,0(\PTR_A,\PTR_B) /*ptrba+off*2**/
			agr \PTR_B, \PTR_B   
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/

		.elseif \C_A==1
			sllg \PTR_B, \OFF_VAL,3 
			agr \PTR_A,\PTR_B /*ptrba+off*4**/
			sllg \PTR_B, \OFF_VAL,5   
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.endif

.elseif \C_B==2
		.if \C_A==8 
			sllg \PTR_B, \OFF_VAL,6 
			agr \PTR_A,\PTR_B /*ptrba+off*8**/
			sllg \PTR_B, \OFF_VAL,4   
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.elseif \C_A==4
			sllg \PTR_B, \OFF_VAL,4 
			la \PTR_A,0(\PTR_A,\PTR_B) /*ptrba+off*2**/
			agr \PTR_A,\PTR_B /*ptrba+off*2**/  
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.elseif \C_A==2
			sllg \PTR_B, \OFF_VAL,4 
			agr \PTR_A,\PTR_B /*ptrba+off*2**/ 
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/ 
		.elseif \C_A==1
			sllg \PTR_B, \OFF_VAL,3 
			la \PTR_A,0(\PTR_A,\PTR_B) /*ptrba+off*1**/
			agr \PTR_B,\PTR_B /* off+off**/  
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.endif

.elseif \C_B==1
		.if \C_A==8 
			sllg \PTR_B, \OFF_VAL,6 
			agr \PTR_A,\PTR_B /*ptrba+off*8**/
			sllg \PTR_B, \OFF_VAL,3 
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.elseif \C_A==4
			sllg \PTR_B, \OFF_VAL,5 
			agr \PTR_A,\PTR_B /*ptrba+off*4**/
			sllg \PTR_B, \OFF_VAL,3 
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.elseif \C_A==2
			sllg \PTR_B, \OFF_VAL,3 
			la \PTR_A,0(\PTR_A,\PTR_B) /*ptrba+off*1**/
			agr \PTR_A,\PTR_B /*ptrba+off*1**/
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/

		.elseif \C_A==1
			sllg \PTR_B, \OFF_VAL,3 
			agr \PTR_A,\PTR_B /*ptrba+off*1**/ 
			la \PTR_B,0(\B_VAL,\PTR_B)	  /*refresh BPOINT*/
		.endif
.endif
             
 
    #endif 
.endm

/**/
.macro RefreshTempBk TEMP_VAL,BK_VAL,OFF_VAL,INCR_A,INCR_B
    #if (defined(LEFT) && !defined(TRANSA)) ||  (!defined(LEFT) && defined(TRANSA))
                            /* temp = bk-off;*/
           sgrk \TEMP_VAL,\BK_VAL,\OFF_VAL

    #elif defined(LEFT)
                            /* temp = off+INCR_A;	// number of values in A */
           la \TEMP_VAL,\INCR_A(\OFF_VAL)
    #else
                            /* temp = off+INCR_B	// number of values in B*/
           la \TEMP_VAL,\INCR_B(\OFF_VAL)
    #endif

.endm
 

.macro RefreshPointersAndOFF TEMP_VAL,BK_VAL,OFF_VAL,PTR_B,PTR_A,C_A,C_B

    #if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
                    /*temp = bk - off;*/
                    sgrk \TEMP_VAL,\BK_VAL,\OFF_VAL
    #ifdef LEFT
                    /*temp -= 8; // number of values in A*/
                    lay \TEMP_VAL,-\C_A(\TEMP_VAL)
    #else
                    /*temp -= 4; // number of values in B*/
                   lay \TEMP_VAL,-\C_B(\TEMP_VAL)
    #endif
                    /*ptrba += temp*C_A;
                    ptrbb += temp*C_B;*/ 
		.if  \C_A==8 
			sllg \TEMP_VAL, \TEMP_VAL,6	
		.elseif \C_A==4
			sllg \TEMP_VAL, \TEMP_VAL,5 /*temp*4*/ 
		.elseif \C_A==2
			sllg \TEMP_VAL, \TEMP_VAL,4 /*temp*2*/  
		.elseif \C_A==1
			sllg \TEMP_VAL, \TEMP_VAL,3 /*temp*1*/  
		.endif
                la \PTR_A,0(\PTR_A,\TEMP_VAL) /*ptrba+temp*C_A*/ 
                /*we do not need to refresh ptrbb. so lets ignore it*/

    #endif

    #ifdef LEFT
                    /*off += 8; // number of values in A*/
                    aghi \OFF_VAL,\C_A
    #endif
.endm