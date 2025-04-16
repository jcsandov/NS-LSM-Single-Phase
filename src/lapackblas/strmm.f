      SUBROUTINE strmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      REAL ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      REAL A(lda,*),B(ldb,*)
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
      REAL TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      REAL ONE,ZERO
      parameter(one=1.0e+0,zero=0.0e+0)
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('STRMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (upper) THEN
                  DO 50 j = 1,n
                      DO 40 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              DO 30 i = 1,k - 1
                                  b(i,j) = b(i,j) + temp*a(i,k)
   30                         CONTINUE
                              IF (nounit) temp = temp*a(k,k)
                              b(k,j) = temp
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              b(k,j) = temp
                              IF (nounit) b(k,j) = b(k,j)*a(k,k)
                              DO 60 i = k + 1,m
                                  b(i,j) = b(i,j) + temp*a(i,k)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A**T*B.
*
              IF (upper) THEN
                  DO 110 j = 1,n
                      DO 100 i = m,1,-1
                          temp = b(i,j)
                          IF (nounit) temp = temp*a(i,i)
                          DO 90 k = 1,i - 1
                              temp = temp + a(k,i)*b(k,j)
   90                     CONTINUE
                          b(i,j) = alpha*temp
  100                 CONTINUE
  110             CONTINUE
              ELSE
                  DO 140 j = 1,n
                      DO 130 i = 1,m
                          temp = b(i,j)
                          IF (nounit) temp = temp*a(i,i)
                          DO 120 k = i + 1,m
                              temp = temp + a(k,i)*b(k,j)
  120                     CONTINUE
                          b(i,j) = alpha*temp
  130                 CONTINUE
  140             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (upper) THEN
                  DO 180 j = n,1,-1
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 150 i = 1,m
                          b(i,j) = temp*b(i,j)
  150                 CONTINUE
                      DO 170 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 160 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  160                         CONTINUE
                          END IF
  170                 CONTINUE
  180             CONTINUE
              ELSE
                  DO 220 j = 1,n
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 190 i = 1,m
                          b(i,j) = temp*b(i,j)
  190                 CONTINUE
                      DO 210 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 200 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
  220             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A**T.
*
              IF (upper) THEN
                  DO 260 k = 1,n
                      DO 240 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              temp = alpha*a(j,k)
                              DO 230 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      temp = alpha
                      IF (nounit) temp = temp*a(k,k)
                      IF (temp.NE.one) THEN
                          DO 250 i = 1,m
                              b(i,k) = temp*b(i,k)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              ELSE
                  DO 300 k = n,1,-1
                      DO 280 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              temp = alpha*a(j,k)
                              DO 270 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  270                         CONTINUE
                          END IF
  280                 CONTINUE
                      temp = alpha
                      IF (nounit) temp = temp*a(k,k)
                      IF (temp.NE.one) THEN
                          DO 290 i = 1,m
                              b(i,k) = temp*b(i,k)
  290                     CONTINUE
                      END IF
  300             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of STRMM .
       END