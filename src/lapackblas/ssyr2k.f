      SUBROUTINE ssyr2k(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      REAL ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      REAL A(lda,*),B(ldb,*),C(ldc,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      REAL TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      REAL ONE,ZERO
      parameter(one=1.0e+0,zero=0.0e+0)
*     ..
*
*     Test the input parameters.
*
      IF (lsame(trans,'N')) THEN
          nrowa = n
      ELSE
          nrowa = k
      END IF
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 1
      ELSE IF ((.NOT.lsame(trans,'N')) .AND.
     +         (.NOT.lsame(trans,'T')) .AND.
     +         (.NOT.lsame(trans,'C'))) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (k.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 7
      ELSE IF (ldb.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldc.LT.max(1,n)) THEN
          info = 12
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('SSYR2K',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR.
     +    (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (upper) THEN
              IF (beta.EQ.zero) THEN
                  DO 20 j = 1,n
                      DO 10 i = 1,j
                          c(i,j) = zero
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 j = 1,n
                      DO 30 i = 1,j
                          c(i,j) = beta*c(i,j)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (beta.EQ.zero) THEN
                  DO 60 j = 1,n
                      DO 50 i = j,n
                          c(i,j) = zero
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 i = j,n
                          c(i,j) = beta*c(i,j)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  C := alpha*A*B**T + alpha*B*A**T + C.
*
          IF (upper) THEN
              DO 130 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 90 i = 1,j
                          c(i,j) = zero
   90                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 100 i = 1,j
                          c(i,j) = beta*c(i,j)
  100                 CONTINUE
                  END IF
                  DO 120 l = 1,k
                      IF ((a(j,l).NE.zero) .OR. (b(j,l).NE.zero)) THEN
                          temp1 = alpha*b(j,l)
                          temp2 = alpha*a(j,l)
                          DO 110 i = 1,j
                              c(i,j) = c(i,j) + a(i,l)*temp1 +
     +                                 b(i,l)*temp2
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 140 i = j,n
                          c(i,j) = zero
  140                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 150 i = j,n
                          c(i,j) = beta*c(i,j)
  150                 CONTINUE
                  END IF
                  DO 170 l = 1,k
                      IF ((a(j,l).NE.zero) .OR. (b(j,l).NE.zero)) THEN
                          temp1 = alpha*b(j,l)
                          temp2 = alpha*a(j,l)
                          DO 160 i = j,n
                              c(i,j) = c(i,j) + a(i,l)*temp1 +
     +                                 b(i,l)*temp2
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*B + alpha*B**T*A + C.
*
          IF (upper) THEN
              DO 210 j = 1,n
                  DO 200 i = 1,j
                      temp1 = zero
                      temp2 = zero
                      DO 190 l = 1,k
                          temp1 = temp1 + a(l,i)*b(l,j)
                          temp2 = temp2 + b(l,i)*a(l,j)
  190                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp1 + alpha*temp2
                      ELSE
                          c(i,j) = beta*c(i,j) + alpha*temp1 +
     +                             alpha*temp2
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 j = 1,n
                  DO 230 i = j,n
                      temp1 = zero
                      temp2 = zero
                      DO 220 l = 1,k
                          temp1 = temp1 + a(l,i)*b(l,j)
                          temp2 = temp2 + b(l,i)*a(l,j)
  220                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp1 + alpha*temp2
                      ELSE
                          c(i,j) = beta*c(i,j) + alpha*temp1 +
     +                             alpha*temp2
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of SSYR2K.
*
      END
