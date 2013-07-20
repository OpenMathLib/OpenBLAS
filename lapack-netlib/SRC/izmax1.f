*> \brief \b IZMAX1 finds the index of the vector element whose real part has maximum absolute value.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download IZMAX1 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/izmax1.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/izmax1.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/izmax1.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       INTEGER          FUNCTION IZMAX1( N, CX, INCX )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         CX( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> IZMAX1 finds the index of the element whose real part has maximum
*> absolute value.
*>
*> Based on IZAMAX from Level 1 BLAS.
*> The change is to use the 'genuine' absolute value.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements in the vector CX.
*> \endverbatim
*>
*> \param[in] CX
*> \verbatim
*>          CX is COMPLEX*16 array, dimension (N)
*>          The vector whose elements will be summed.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The spacing between successive values of CX.  INCX >= 1.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*> Nick Higham for use with ZLACON.
*
*  =====================================================================
      INTEGER          FUNCTION IZMAX1( N, CX, INCX )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX
      DOUBLE PRECISION   SMAX
      COMPLEX*16         ZDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
*
*     NEXT LINE IS THE ONLY MODIFICATION.
      CABS1( ZDUM ) = ABS( ZDUM )
*     ..
*     .. Executable Statements ..
*
      IZMAX1 = 0
      IF( N.LT.1 )
     $   RETURN
      IZMAX1 = 1
      IF( N.EQ.1 )
     $   RETURN
      IF( INCX.EQ.1 )
     $   GO TO 30
*
*     CODE FOR INCREMENT NOT EQUAL TO 1
*
      IX = 1
      SMAX = CABS1( CX( 1 ) )
      IX = IX + INCX
      DO 20 I = 2, N
         IF( CABS1( CX( IX ) ).LE.SMAX )
     $      GO TO 10
         IZMAX1 = I
         SMAX = CABS1( CX( IX ) )
   10    CONTINUE
         IX = IX + INCX
   20 CONTINUE
      RETURN
*
*     CODE FOR INCREMENT EQUAL TO 1
*
   30 CONTINUE
      SMAX = CABS1( CX( 1 ) )
      DO 40 I = 2, N
         IF( CABS1( CX( I ) ).LE.SMAX )
     $      GO TO 40
         IZMAX1 = I
         SMAX = CABS1( CX( I ) )
   40 CONTINUE
      RETURN
*
*     End of IZMAX1
*
      END
