      DOUBLE PRECISION FUNCTION NEWTON_P(Q,LDQ,T,LDT,SIG,LDSIG,Y,
     $                                   LDY,N,S,ITERATION)
      IMPLICIT NONE
!variables données :
      INTEGER LDQ,LDT,LDSIG,LDY,N,ITERATION,I,J
      DOUBLE PRECISION Q(LDQ,*),T(LDT,*),Y(LDY),SIG(LDSIG,*),S
!variables locales :
      DOUBLE PRECISION QTQ(N-1,N-1),M(N-1,N-1),INVERSE(N-1,N-1)
      DOUBLE PRECISION QQT(N-1,N+1),SIGR(N+1,N+1),U(N-1),V(N+1)
      DOUBLE PRECISION SIGQ(N+1,N-1)
      DOUBLE PRECISION UMAT(1,N-1), QQ(1,N-1),QSIG(N-1,N+1)
      DOUBLE PRECISION FP2, FP2P, NORME2, FP, FPP
      NEWTON_P = 0.0001D0
      DO I=1,N+1
        DO J=1,N+1
            IF (I==J) THEN
                SIGR(I,I) = SQRT(SIG(I,I))
            ELSE
                SIGR(I,J) = 0
            END IF
        END DO
      END DO
      DO I=1,ITERATION
        CALL DGEMM('T','N',N-1,LDSIG,LDQ,1D0,Q,LDQ,SIG,LDSIG,
     $           0D0,QSIG,N-1)
        CALL DGEMM('N','N',N-1,N-1,N+1,1D0,QSIG,N-1,Q,LDQ,
     $           0D0,QTQ,N-1)
        CALL ADDITION_MATRICE_CARRE(QTQ,N-1,T,LDT,M,N-1,NEWTON_P)
        CALL INVERSION(M,N-1,INVERSE,N-1,N-1)
        CALL DGEMM('N','T',N-1,N+1,N-1,1D0,INVERSE,N-1,Q,LDQ,
     $           0D0,QQT,N-1)
        CALL DGEMV('N',N-1,N+1,1D0,QQT,N-1,Y,1,0D0,U,1)
        CALL DGEMM('N','N',N+1,N-1,N+1,1D0,SIGR,N+1,Q,LDQ,0D0,
     $             SIGQ,N+1)
        CALL DGEMV('N',N+1,N-1,1D0,SIGQ,N+1,U,1,0D0,V,1)
        FP2 = NORME2(V,N+1,V,N+1)
        UMAT(1,:) = U
        CALL DGEMM('N','N',1,N-1,N-1,1D0,UMAT,1,T,LDT,0D0,
     $             QQ,1)
        UMAT(1,:) = QQ(1,:)
        FP2P = -1D0 * NORME2(UMAT,N-1,U,N-1)
        CALL DGEMM('N','N',1,N-1,N-1,1D0,UMAT,1,INVERSE,N-1,0D0,
     $             QQ,1)
        UMAT(1,:) = QQ(1,:)
        CALL DGEMM('N','N',1,N-1,N-1,1D0,UMAT,1,T,LDT,0D0,
     $             QQ,1)
        FP2P = 2D0 * (NEWTON_P * NORME2(QQ,N-1,U,N-1) + FP2P)
        FP = (1/SQRT(FP2)) - (1/SQRT(S))
        FPP = (-0.5D0 * FP2P * (1/SQRT(FP2)))/FP2
        NEWTON_P = NEWTON_P - (FP/FPP)
      END DO
      WRITE(*,*) "p -> ",NEWTON_P
      RETURN
      END FUNCTION
