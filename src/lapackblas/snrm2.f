*> \brief \b SNRM2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       REAL FUNCTION SNRM2(N,X,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       REAL X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    SNRM2 := sqrt( x'*x ).
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
*> \date December 2016
*
*> \ingroup single_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  -- This version written on 25-October-1982.
*>     Modified on 14-October-1993 to inline the call to SLASSQ.
*>     Sven Hammarling, Nag Ltd.
*> \endverbatim
*>
*  =====================================================================
      REAL FUNCTION snrm2(N,X,INCX)
*
*  -- Reference BLAS level1 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL ONE,ZERO
      parameter(one=1.0e+0,zero=0.0e+0)
*     ..
*     .. Local Scalars ..
      REAL ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt
*     ..
      IF (n.LT.1 .OR. incx.LT.1) THEN
          norm = zero
      ELSE IF (n.EQ.1) THEN
          norm = abs(x(1))
      ELSE
          scale = zero
          ssq = one
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 ix = 1,1 + (n-1)*incx,incx
              IF (x(ix).NE.zero) THEN
                  absxi = abs(x(ix))
                  IF (scale.LT.absxi) THEN
                      ssq = one + ssq* (scale/absxi)**2
                      scale = absxi
                  ELSE
                      ssq = ssq + (absxi/scale)**2
                  END IF
              END IF
   10     CONTINUE
          norm = scale*sqrt(ssq)
      END IF
*
      snrm2 = norm
      RETURN
*
*     End of SNRM2.
*
      END
